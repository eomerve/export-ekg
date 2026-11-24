#!/usr/bin/env python3
import argparse, csv, json, struct, os
from pathlib import Path
from bisect import bisect_left
os.umask(0o077)
p=argparse.ArgumentParser(description='BT740 research-only waveform candidate export')
p.add_argument('--capture', type=Path, default=Path.home()/'bistos-capture/bistos-capture/capture_20261124_154913')
p.add_argument('--offset',type=int)
p.add_argument('--count',type=int)
p.add_argument('--fs',type=float)
p.add_argument('--stride',type=int,default=2)
p.add_argument('--resp-offset',type=int)
p.add_argument('--resp-count',type=int)
p.add_argument('--resp-stride',type=int,default=2)
p.add_argument('--resp-fs',type=float)
p.add_argument('--gain',type=float)
p.add_argument('--baseline',type=int,default=0)
a=p.parse_args()
raw=(a.capture/'cms_raw.bin').read_bytes()
chunks=[json.loads(l) for l in (a.capture/'cms_raw.txt').read_text().splitlines() if l.strip()]
ends=[c['byte_offset']+c['length'] for c in chunks]
frames=[]; pos=0
while pos+10<=len(raw):
    if raw[pos:pos+4]!=b'QCMS': raise SystemExit(f'Invalid header at {pos}')
    n=int.from_bytes(raw[pos+8:pos+10],'little')
    if n<10 or pos+n>len(raw): print('Incomplete frame at',pos); break
    b=raw[pos:pos+n]
    if n==2518 and b[7]==5:
        i=bisect_left(ends,pos+n)
        t=chunks[i]['timestamp'] if i<len(chunks) and chunks[i]['byte_offset']<pos+n else None
        frames.append((b,t))
    pos+=n
if not frames: raise SystemExit('No periodic frames')
out=a.capture/'ecg_analysis'; out.mkdir(exist_ok=True)
print('Frames:',len(frames))
for i,(_,t) in enumerate(frames): print(i,t)
# Candidate blocks, not identified signals. Both byte alignments.
for alignment in (0,1):
    with (out/f'candidate_{alignment}.csv').open('w',newline='') as f:
        w=csv.writer(f); w.writerow(['frame','timestamp','byte_offset','raw_int16'])
        for i,(b,t) in enumerate(frames):
            for j in range(64+alignment,len(b)-1,2): w.writerow([i,t,j,struct.unpack_from('<h',b,j)[0]])
# List changing blocks with value range, to narrow waveform search without plots.
for s in range(64,2518,64):
    e=min(s+64,2518)
    v=[struct.unpack_from('<h',b,j)[0] for b,_ in frames for j in range(s,e-1,2)]
    changed=sum(len({b[j] for b,_ in frames})>1 for j in range(s,e))
    print(f'{s:4}-{e-1:4} changed={changed:2} min={min(v):6} max={max(v):6}')
if a.offset is None:
    print('Candidate CSV saved. No verified ECG export yet.'); raise SystemExit(0)
def validate(o,c,s,fs):
    if o is None or c is None or fs is None or o<16 or c<1 or s<2 or fs<=0 or o+(c-1)*s+2>2518: raise SystemExit('Invalid or missing signal layout/rate')
validate(a.offset,a.count,a.stride,a.fs)
resp=any(x is not None for x in (a.resp_offset,a.resp_count,a.resp_fs))
if resp: validate(a.resp_offset,a.resp_count,a.resp_stride,a.resp_fs)
if a.gain is not None and a.gain<=0: raise SystemExit('Gain must be positive')
# Keep experimental outputs separate; never append across captures.
out=a.capture/'experimental_export'; out.mkdir(exist_ok=False)
def samples(b,o,c,s): return [struct.unpack_from('<h',b,o+i*s)[0] for i in range(c)]
all_ecg=[]; all_resp=[]
with (out/'data.txt').open('w') as f:
    for b,t in frames:
        x=samples(b,a.offset,a.count,a.stride)
        d={'timestamp':t,'status':'experimental_unverified','vitals_candidates':{'HR_offset56':struct.unpack_from('<h',b,56)[0],'PR_offset26':struct.unpack_from('<h',b,26)[0],'PR_offset36':struct.unpack_from('<h',b,36)[0]},'ecg':{'start_sample':len(all_ecg),'raw':x}}
        all_ecg.extend(x)
        if resp:
            r=samples(b,a.resp_offset,a.resp_count,a.resp_stride); d['resp']={'start_sample':len(all_resp),'raw':r}; all_resp.extend(r)
        f.write(json.dumps(d)+'\n')
(out/'ecg.dat').write_bytes(struct.pack('<'+'h'*len(all_ecg),*all_ecg))
gain=a.gain or 1; unit='mV' if a.gain is not None else 'adu'
checksum=((sum(all_ecg)+32768)%65536)-32768
(out/'ecg.hea').write_text(f'ecg 1 {a.fs:g} {len(all_ecg)}\necg.dat 16 {gain:g}({a.baseline})/{unit} 16 0 {all_ecg[0]} {checksum} 0 ECG_candidate\n# Experimental unverified layout and sample rate; not clinical data\n')
(out/'ecg.atr').write_bytes(b'\x00\x00')
if resp: (out/'ecg.resp').write_bytes(struct.pack('<'+'h'*len(all_resp),*all_resp))
(out/'metadata.json').write_text(json.dumps({'experimental':True,'ecg_fs_hz':a.fs,'ecg_offset':a.offset,'ecg_stride':a.stride,'ecg_count_per_frame':a.count,'resp_fs_hz':a.resp_fs,'resp_offset':a.resp_offset,'resp_stride':a.resp_stride,'resp_count_per_frame':a.resp_count,'resp_encoding':'signed int16 little endian' if resp else None,'timing':'TCP receive completion timestamps; gaps not reconstructed','calibration_verified':False},indent=2))
print('Experimental export:',out)
