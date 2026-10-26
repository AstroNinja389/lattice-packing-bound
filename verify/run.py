"""Run one independent verifier and compare its results with reference logs."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
VERIFIERS=('pari','taylor','bernstein','descartes')


def normalized(output):
    return '\n'.join(line.rstrip() for line in output.splitlines()
                     if not line.lower().startswith(('time:','seconds:','elapsed:'))).strip()


def unique_object(pairs):
    result={}
    for key,value in pairs:
        if key in result:raise ValueError('duplicate JSON key')
        result[key]=value
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--verifier',choices=VERIFIERS,required=True)
    p.add_argument('--file',nargs='+',type=Path)
    p.add_argument('--manifest',type=Path,default=ROOT/'manifest.json')
    p.add_argument('--gp',default=os.environ.get('PARI_GP'))
    p.add_argument('--timeout',type=float,default=3600)
    p.add_argument('--jobs',type=int,default=1,help='independent verifier processes, 1..16 (default: 1)')
    p.add_argument('--no-reference',action='store_true',help='verify modified copies without comparing reference output')
    p.add_argument('--full-run',action='store_true',help='explicitly authorize recording the complete manifest run')
    p.add_argument('--record-reference',action='store_true')
    a=p.parse_args()
    if not math.isfinite(a.timeout) or a.timeout<=0:p.error('timeout must be finite and positive')
    if not 1<=a.jobs<=16:p.error('jobs must lie between 1 and 16')
    if a.record_reference and (not a.full_run or a.file or a.no_reference):p.error('--record-reference requires --full-run and the complete manifest')
    if a.full_run and a.file:p.error('--full-run cannot be combined with --file')
    references=ROOT/'logs'/f'{a.verifier}.json'
    try:
        if a.file:
            paths=a.file
        else:
            manifest=json.loads(a.manifest.read_text(encoding='utf-8'),object_pairs_hook=unique_object)
            entries=manifest['certificates']
            if type(entries) is not list or not entries:raise ValueError('empty manifest')
            paths=[];names=set()
            for item in entries:
                path=(a.manifest.parent/item['file']).resolve()
                if not path.is_relative_to(a.manifest.parent.resolve()):raise ValueError('manifest path escapes repository')
                if path.name.casefold() in names:raise ValueError('duplicate manifest certificate')
                names.add(path.name.casefold())
                if hashlib.sha256(path.read_bytes()).hexdigest()!=item['sha256']:raise ValueError('manifest hash mismatch: '+path.name)
                header=json.loads(path.read_text(encoding='utf-8'),object_pairs_hook=unique_object)
                for key in ('n','K'):
                    if type(header.get(key)) is not int or type(item.get(key)) is not int or header[key]!=item[key]:
                        raise ValueError('manifest metadata mismatch: '+path.name)
                if item.get('claims') != header.get('claims'):
                    raise ValueError('manifest claims mismatch: '+path.name)
                if a.verifier=='descartes' and header['n']%2:continue
                paths.append(path)
        if not paths:raise ValueError('no applicable certificates to check')
        if len({path.name.casefold() for path in paths})!=len(paths):raise ValueError('duplicate certificate basename')
        reference={}
        if not a.no_reference and not a.record_reference:
            reference=json.loads(references.read_text(encoding='utf-8'))['results']
        def run_one(path):
            command=[sys.executable,str(ROOT/'verify'/f'{a.verifier}.py'),str(path.resolve())]
            if a.verifier=='pari' and a.gp:command+=['--gp',a.gp]
            start=time.perf_counter()
            try:
                run=subprocess.run(command,capture_output=True,text=True,timeout=a.timeout)
                return path,run,round(time.perf_counter()-start,6),None
            except (OSError,subprocess.SubprocessError) as exc:
                return path,None,round(time.perf_counter()-start,6),str(exc)
        results={};timings={};passed=True
        with ThreadPoolExecutor(max_workers=a.jobs) as pool:
            for path,run,elapsed,error in pool.map(run_one,paths):
                if error is not None:
                    passed=False;print(path.name+': '+error,file=sys.stderr);continue
                output=normalized(run.stdout);timings[path.name]=elapsed
                good=run.returncode==0 and output.startswith('PASS ')
                if not a.no_reference and not a.record_reference:
                    good=good and reference.get(path.name)==output
                if good:
                    results[path.name]=output;print(path.name+': '+output,flush=True)
                else:
                    passed=False;print(path.name+': REJECT or reference mismatch',file=sys.stderr)
                    print((run.stdout+run.stderr).strip(),file=sys.stderr)
        if not passed or len(results)!=len(paths):return 1
        if a.record_reference:
            references.parent.mkdir(parents=True,exist_ok=True)
            payload={'verifier':a.verifier,'results':results,'timings_seconds':timings}
            temporary=references.with_suffix('.tmp');temporary.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n',encoding='utf-8');temporary.replace(references)
        print(f'[computed] Checked {len(results)} certificate(s); all passed.')
        return 0
    except (OSError,ValueError,KeyError,TypeError) as exc:
        print('[computed] REJECT: '+str(exc),file=sys.stderr);return 1


if __name__=='__main__':sys.exit(main())
