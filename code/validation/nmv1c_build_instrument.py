"""NM-V1C -- fresh draw and Method-X-free offline adjudication application.

Method X is not merely hidden: it is never loaded. mge_feature_inventory.tsv, the HMM feature
source, is not opened by this script at all, so no Method X coordinate can reach the instrument
through any path including hidden markup, script state, metadata or exports.
"""
import argparse, collections, csv, datetime, hashlib, html, json, os, sys

VERSION = "nmv1c_build_instrument_v1.1.0"
PROTO_SHA = "b2058877c0f7165f0ab7bf9be7adcb2cb3b5b4c8d2b81143d85ab03fe0f0c04c"
REASONS = ["INCOMPLETE_ELEMENT", "AMBIGUOUS_EVIDENCE", "BOUNDARY_TRUNCATED",
           "TOOL_FAILURE", "INSUFFICIENT_EVIDENCE"]


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def dsort(x, seed, salt):
    return sorted(x, key=lambda b: hashlib.sha256(
        ("%s|%d|%s" % (b, seed, salt)).encode()).hexdigest())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--dir", required=True)
    ap.add_argument("--pool", required=True)
    ap.add_argument("--suffix", default="",
                    help="written to *<suffix>.html/.json/.tsv; never overwrites an "
                         "existing build")
    a = ap.parse_args()
    SFX = a.suffix
    D = a.dir
    proto = os.path.join(D, "NMV1C_FROZEN_PROTOCOL.json")
    if sha(proto) != PROTO_SHA:
        print("REFUSING: protocol digest mismatch"); sys.exit(1)
    for f in ("NMV1C_ADJUDICATION_APP%s.html" % SFX,
              "NMV1C_SAMPLE_MANIFEST%s.json" % SFX,
              "NMV1C_UNBLINDING_KEY%s.tsv" % SFX):
        if os.path.exists(os.path.join(D, f)):
            print("REFUSING: %s already exists; a build never overwrites one" % f); sys.exit(1)
    P = json.load(open(proto, encoding="utf-8"))
    if not P.get("frozen_before_any_case_was_drawn"):
        print("REFUSING: protocol does not assert pre-draw freeze"); sys.exit(1)
    SEED = P["selection"]["seed"]; ALLOC = P["selection"]["allocation"]
    POOL = json.load(open(a.pool, encoding="utf-8"))
    excluded = set(POOL["excluded"]); avail = POOL["available"]
    print("%s | protocol %s verified" % (VERSION, PROTO_SHA[:16]))

    gt = {r["block_id"]: r for r in csv.DictReader(
        open(os.path.join(D, "NMV1_RULE_BASED_GROUND_TRUTH.tsv"), encoding="utf-8"),
        delimiter="\t")}
    ev = {r["block_id"]: r for r in csv.DictReader(
        open(os.path.join(D, "nmv1_block_evidence_table.tsv"), encoding="utf-8"),
        delimiter="\t")}

    # ---------------- fresh stratified draw ----------------
    sel = {}; used = set()
    for st in ("A", "B", "D", "E"):
        cand = sorted([b for b in avail if gt[b]["rule_id"] == st and b not in used])
        pick = dsort(cand, SEED, st)[:ALLOC[st]]
        sel[st] = pick; used.update(pick)
        print("  stratum %-2s pool %4d  allocated %3d  drawn %3d  weight %.4f"
              % (st, len(cand), ALLOC[st], len(pick), len(cand) / max(len(pick), 1)))
    fresh = sorted(used)
    assert len(fresh) == 120, len(fresh)
    assert not (set(fresh) & excluded), "overlap with a previously shown block"
    print("  fresh sample: %d blocks | overlap with any prior package: 0" % len(fresh))

    weight = {b: len([x for x in avail if gt[x]["rule_id"] == gt[b]["rule_id"]]) / ALLOC[gt[b]["rule_id"]]
              for b in fresh}

    # ---------------- evidence, Method Y and Z ONLY ----------------
    O = os.path.join(a.repo, "audit/data/derived/pr_context/out")
    RES = os.path.join(a.repo, "audit/data/derived/nmv1")
    blk = {r["block_id"]: r for r in csv.DictReader(
        open(os.path.join(O, "shared_context_blocks.tsv"), encoding="utf-8"), delimiter="\t")
        if r["block_id"] in set(fresh)}
    ise = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(RES, "isescan_hits.tsv"), encoding="utf-8"),
                            delimiter="\t"):
        if r["block_id"] in blk:
            ise[r["block_id"]].append(r)
    inf = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(RES, "integronfinder_features.tsv"),
                                 encoding="utf-8"), delimiter="\t"):
        if r["block_id"] in blk:
            inf[r["block_id"]].append(r)
    # NOTE: mge_feature_inventory.tsv is deliberately never opened.
    tgt = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(O, "arg_mge_neighbourhood.tsv"),
                                 encoding="utf-8"), delimiter="\t"):
        for b, s in blk.items():
            if s["replicon_accession"] != r["replicon_accession"]:
                continue
            bs, be = int(s["block_start"]), int(s["block_end"])
            g1, g2 = int(r["gene_start"]), int(r["gene_end"])
            if bs <= g1 and be >= g2:
                tgt[b].append((g1 - bs + 1, g2 - bs + 1))
                break

    def num(x):
        try:
            return int(float(x))
        except (TypeError, ValueError):
            return 0

    tok = {b: "NMC-%03d" % (i + 1) for i, b in enumerate(dsort(fresh, SEED, "nmv1ctoken"))}
    cases = []
    for b in dsort(fresh, SEED, "nmv1ctoken"):
        s = blk[b]; L = int(s["block_span_bp"]); T = tgt.get(b, [])
        F = []
        for h in ise.get(b, []):
            ir = ""
            if num(h.get("irLen")) > 0:
                ir = "left %s-%s / right %s-%s, %s bp, identity %s" % (
                    h["start1"], h["end1"], h["start2"], h["end2"], h["irLen"], h["irId"])
            F.append({"m": "Y", "t": "complete element" if h.get("type") == "c" else "partial element",
                      "b": num(h["isBegin"]), "e": num(h["isEnd"]), "s": h.get("strand", ""),
                      "v": h.get("E_value", ""),
                      "d": "inverted repeats: %s | ORF %s-%s, %s bp"
                           % (ir or "not resolved", h.get("orfBegin", ""), h.get("orfEnd", ""),
                              h.get("orfLen", ""))})
        for f in inf.get(b, []):
            pb, pe = num(f["pos_beg"]), num(f["pos_end"])
            if not pb and not pe:
                continue
            t = ("attC site" if f.get("type_elt") == "attC"
                 else "integrase" if f.get("annotation") == "intI" else "cassette protein")
            F.append({"m": "Z", "t": t, "b": pb, "e": pe, "s": f.get("strand", ""),
                      "v": f.get("evalue", ""),
                      "d": "structure: %s" % f.get("type", "")})
        F.sort(key=lambda z: z["b"])
        cases.append({"token": tok[b], "len": L, "topology": s["topology"],
                      "boundary": "yes" if (s["truncated"] == "yes"
                                            or s["wrapped_circular"] == "yes") else "no",
                      "tool": ev[b]["tool_status"], "target": T, "features": F})
    ny = sum(1 for c in cases if any(f["m"] == "Y" for f in c["features"]))
    nz = sum(1 for c in cases if any(f["m"] == "Z" for f in c["features"]))
    print("  cases with Method Y evidence: %d | Method Z: %d | no evidence: %d"
          % (ny, nz, sum(1 for c in cases if not c["features"])))

    # ---------------- offline application ----------------
    app = os.path.join(D, "NMV1C_ADJUDICATION_APP%s.html" % SFX)
    payload = json.dumps(cases, separators=(",", ":"))
    H = []
    H.append("<!doctype html><html lang='en'><meta charset='utf-8'>")
    H.append("<title>NM-V1C blinded adjudication</title>")
    H.append("<style>"
             ":root{--bg:#fff;--fg:#1a1a1a;--mut:#666;--line:#ddd;--warn:#8a4b00}"
             "*{box-sizing:border-box}body{margin:0;font:15px/1.55 ui-sans-serif,system-ui,sans-serif;"
             "background:var(--bg);color:var(--fg)}"
             "header{position:sticky;top:0;background:#fff;border-bottom:1px solid var(--line);"
             "padding:10px 20px;display:flex;gap:18px;align-items:center;flex-wrap:wrap;z-index:9}"
             "#bar{flex:1;height:8px;background:#eee;border-radius:4px;overflow:hidden;min-width:160px}"
             "#fill{height:100%;background:#5d8aa8;width:0}"
             "main{max-width:1040px;margin:0 auto;padding:20px}"
             "h1{font-size:17px;margin:0}"
             ".meta{color:var(--mut);font-size:13px;margin:6px 0 14px}"
             ".warn{color:var(--warn);font-weight:600}"
             "table{border-collapse:collapse;font-size:12.5px;width:100%;margin:10px 0}"
             "th,td{border:1px solid var(--line);padding:4px 8px;text-align:left}"
             "th{background:#fafafa}"
             ".choices{display:flex;gap:10px;margin:14px 0 6px;flex-wrap:wrap}"
             "button.c{padding:10px 16px;border:2px solid var(--line);background:#fff;"
             "border-radius:8px;cursor:pointer;font:inherit}"
             "button.c[aria-pressed=true]{border-color:#1b2631;background:#eef3f7;font-weight:700}"
             "kbd{border:1px solid #bbb;border-bottom-width:2px;border-radius:4px;padding:0 5px;"
             "font-size:11px;color:#444}"
             "select,textarea{font:inherit;padding:7px;border:1px solid var(--line);"
             "border-radius:6px;width:100%;max-width:640px}"
             "nav{display:flex;gap:10px;margin:22px 0}"
             "nav button{padding:8px 14px;border:1px solid var(--line);background:#fff;"
             "border-radius:6px;cursor:pointer;font:inherit}"
             "#req{color:#b00;font-size:13px;display:none;margin-top:6px}"
             ".legend span{display:inline-block;width:11px;height:11px;border:1px solid #444;"
             "vertical-align:-1px;margin-right:4px}"
             ".legend{font-size:12px;color:var(--mut);margin-top:8px}"
             ".done{color:#2a6b2a;font-weight:600}"
             ".imp{font-size:12px;color:var(--mut);cursor:pointer}"
             ".imp input{display:block;font-size:11px;max-width:190px}"
             "#nostore{display:none;margin:0;padding:12px 20px;background:#fff4e5;"
             "border-bottom:1px solid #e0c9a6;color:var(--warn);font-size:13.5px}"
             "</style>")
    H.append("<header><h1>NM-V1C blinded adjudication</h1>"
             "<div id='bar'><div id='fill'></div></div>"
             "<span id='prog' class='meta' style='margin:0'></span>"
             "<button onclick='exp(\"tsv\")'>Export TSV</button>"
             "<button onclick='exp(\"json\")'>Export JSON</button>"
             "<label class='imp'>Restore progress"
             "<input type='file' accept='.json,application/json' "
             "onchange='if(this.files[0])imp(this.files[0])'></label>"
             "</header>")
    H.append("<div id='nostore' role='alert'></div>")
    H.append("<main>")
    H.append("<p class='meta'>Judge only the structural evidence shown. Two independent "
             "detection methods are shown as <b>Method Y</b> and <b>Method Z</b>. Species, "
             "study, gene identity and every machine-derived label are withheld. Coordinates "
             "are 1-based and relative to the block. Decisions autosave in this browser.</p>")
    H.append("<p class='meta'><b>MOBILE</b> — independent structural evidence supports a "
             "mobile-element context. &nbsp; <b>QUIESCENT</b> — independent tools completed "
             "normally and no credible mobile-element evidence is present. &nbsp; "
             "<b>NON_EVALUABLE</b> — independent evidence is incomplete, ambiguous, "
             "boundary-truncated, technically failed or insufficient to establish or exclude "
             "mobile context. A reason code is required for NON_EVALUABLE.</p>")
    H.append("<p class='meta'>Keyboard: <kbd>1</kbd> mobile &nbsp; <kbd>2</kbd> quiescent "
             "&nbsp; <kbd>3</kbd> non-evaluable &nbsp; <kbd>&larr;</kbd>/<kbd>&rarr;</kbd> or "
             "<kbd>j</kbd>/<kbd>k</kbd> navigate</p>")
    H.append("<div id='case'></div>")
    H.append("<nav><button onclick='go(-1)'>&larr; Previous</button>"
             "<button onclick='go(1)'>Next &rarr;</button>"
             "<button onclick='nextUnanswered()'>Next unanswered</button></nav>")
    H.append("<p class='legend'><span style='background:#9e9e9e'></span>target interval "
             "<span style='background:#5d8aa8'></span>Method Y complete element "
             "<span style='background:#a9cce3'></span>Method Y partial element "
             "<span style='background:#f0b27a'></span>Method Z attC "
             "<span style='background:#c39bd3'></span>Method Z integrase "
             "<span style='background:#d5dbdb'></span>Method Z cassette protein</p>")
    H.append("</main>")
    H.append("<script>const CASES=%s;const REASONS=%s;" % (payload, json.dumps(REASONS)))
    H.append(r"""
const KEY="nmv1c_decisions_v1";
let STORE_OK=true, STORE_WHY="";
function lget(){try{return localStorage.getItem(KEY);}
  catch(e){STORE_OK=false;STORE_WHY=e.message||String(e);return null;}}
function lset(v){try{localStorage.setItem(KEY,v);return true;}
  catch(e){if(STORE_OK){STORE_OK=false;STORE_WHY=e.message||String(e);banner();}return false;}}
let i=0, dec={};
try{dec=JSON.parse(lget()||"{}")||{};}catch(e){dec={};}
const col={"complete element":"#5d8aa8","partial element":"#a9cce3","attC site":"#f0b27a",
           "integrase":"#c39bd3","cassette protein":"#d5dbdb"};
function svg(c){
  const W=960,H=150,sc=x=>40+(W-80)*Math.max(0,Math.min(c.len,x))/Math.max(c.len,1);
  let p=`<svg width="100%" viewBox="0 0 ${W} ${H}" font-family="ui-monospace,monospace" font-size="11">`;
  p+=`<rect width="${W}" height="${H}" fill="#fff"/>`;
  p+=`<line x1="${sc(0)}" y1="26" x2="${sc(c.len)}" y2="26" stroke="#333"/>`;
  for(let k=0;k<=5;k++){const x=sc(c.len*k/5);
    p+=`<line x1="${x}" y1="22" x2="${x}" y2="30" stroke="#333"/>`;
    p+=`<text x="${x}" y="17" text-anchor="middle" fill="#333">${Math.round(c.len*k/5)}</text>`;}
  p+=`<text x="6" y="17" fill="#333">bp</text>`;
  const tr={"T":46,"Y":76,"Z":106};
  for(const [k,y] of Object.entries({"target":46,"Method Y":76,"Method Z":106})){
    p+=`<text x="6" y="${y+10}" fill="#555">${k}</text>`;
    p+=`<line x1="${sc(0)}" y1="${y+7}" x2="${sc(c.len)}" y2="${y+7}" stroke="#eee"/>`;}
  (c.target||[]).forEach(t=>{p+=`<rect x="${sc(t[0])}" y="46" width="${Math.max(2,sc(t[1])-sc(t[0]))}" height="12" fill="#9e9e9e" stroke="#555"/>`;});
  (c.features||[]).forEach(f=>{const y=f.m=="Y"?76:106;const x1=sc(f.b),x2=sc(f.e);
    p+=`<rect x="${x1}" y="${y}" width="${Math.max(2,x2-x1)}" height="12" fill="${col[f.t]||"#ccc"}" stroke="#444"><title>${f.t} ${f.b}-${f.e}</title></rect>`;
    if(f.t=="complete element"&&/left \d/.test(f.d)){[x1,x2].forEach(x=>{p+=`<rect x="${x-2}" y="${y}" width="4" height="12" fill="#1b2631"/>`;});}});
  return p+"</svg>";}
function esc(s){return String(s).replace(/[&<>]/g,m=>({"&":"&amp;","<":"&lt;",">":"&gt;"}[m]));}
function render(){
  const c=CASES[i], d=dec[c.token]||{};
  let h=`<h2 style="font-size:16px;margin:0 0 4px">${c.token} <span class="meta" style="font-weight:400">(${i+1} of ${CASES.length})</span></h2>`;
  let w="";
  if(c.boundary=="yes") w+=` <span class="warn">boundary warning: element may extend beyond the retrieved sequence</span>`;
  if(c.tool!="ok") w+=` <span class="warn">independent tool status: ${esc(c.tool)}</span>`;
  h+=`<div class="meta">length <b>${c.len} bp</b> &middot; topology <b>${esc(c.topology)}</b> &middot; target interval(s) <b>${(c.target||[]).map(t=>t[0]+"-"+t[1]).join("; ")||"none"}</b>${w}</div>`;
  h+=svg(c);
  if(c.features.length){
    h+=`<table><tr><th>method</th><th>feature</th><th>start</th><th>end</th><th>strand</th><th>e-value</th><th>structural detail</th></tr>`;
    c.features.forEach(f=>{h+=`<tr><td>${f.m}</td><td>${esc(f.t)}</td><td>${f.b}</td><td>${f.e}</td><td>${esc(f.s)}</td><td>${esc(f.v)}</td><td>${esc(f.d)}</td></tr>`;});
    h+=`</table>`;
  } else { h+=`<p class="meta"><i>Neither independent method reported any feature in this block.</i></p>`; }
  h+=`<div class="choices">
    <button class="c" id="bM" aria-pressed="${d.decision=="MOBILE"}" onclick="pick('MOBILE')">MOBILE <kbd>1</kbd></button>
    <button class="c" id="bQ" aria-pressed="${d.decision=="QUIESCENT"}" onclick="pick('QUIESCENT')">QUIESCENT <kbd>2</kbd></button>
    <button class="c" id="bN" aria-pressed="${d.decision=="NON_EVALUABLE"}" onclick="pick('NON_EVALUABLE')">NON_EVALUABLE <kbd>3</kbd></button></div>`;
  h+=`<div id="rwrap" style="display:${d.decision=="NON_EVALUABLE"?"block":"none"}">
      <label class="meta">Reason code (required for NON_EVALUABLE)</label>
      <select id="rc" onchange="setr(this.value)"><option value="">-- select --</option>
      ${REASONS.map(r=>`<option value="${r}" ${d.reason==r?"selected":""}>${r}</option>`).join("")}</select>
      <div id="req">A reason code is required before this case counts as answered.</div></div>`;
  h+=`<label class="meta" style="display:block;margin-top:12px">Note (optional)</label>
      <textarea id="nt" rows="2" oninput="setn(this.value)">${esc(d.note||"")}</textarea>`;
  document.getElementById("case").innerHTML=h;
  const done=CASES.filter(x=>ok(dec[x.token])).length;
  document.getElementById("prog").innerHTML=`<span class="${done==CASES.length?"done":""}">${done} / ${CASES.length} answered</span>`;
  document.getElementById("fill").style.width=(100*done/CASES.length)+"%";
  if(d.decision=="NON_EVALUABLE"&&!d.reason) document.getElementById("req").style.display="block";
}
function ok(d){return d&&d.decision&&(d.decision!="NON_EVALUABLE"||d.reason);}
function save(){lset(JSON.stringify(dec));}
function pick(v){const t=CASES[i].token;dec[t]=Object.assign({},dec[t],{decision:v});
  if(v!="NON_EVALUABLE")delete dec[t].reason;save();render();}
function setr(v){const t=CASES[i].token;dec[t]=Object.assign({},dec[t],{reason:v});save();render();}
function setn(v){const t=CASES[i].token;dec[t]=Object.assign({},dec[t],{note:v});save();}
function go(n){i=Math.max(0,Math.min(CASES.length-1,i+n));render();window.scrollTo(0,0);}
function nextUnanswered(){for(let k=1;k<=CASES.length;k++){const j=(i+k)%CASES.length;
  if(!ok(dec[CASES[j].token])){i=j;render();window.scrollTo(0,0);return;}}alert("All cases answered.");}
document.addEventListener("keydown",e=>{
  if(["INPUT","TEXTAREA","SELECT"].includes(e.target.tagName))return;
  if(e.key=="1")pick("MOBILE");else if(e.key=="2")pick("QUIESCENT");
  else if(e.key=="3")pick("NON_EVALUABLE");
  else if(e.key=="ArrowRight"||e.key=="j")go(1);else if(e.key=="ArrowLeft"||e.key=="k")go(-1);});
function exp(fmt){
  const miss=CASES.filter(c=>!ok(dec[c.token]));
  if(miss.length&&!confirm(miss.length+" case(s) are unanswered or missing a reason code. Export anyway?"))return;
  let blob,name;
  if(fmt=="tsv"){let s="token\tdecision\treason_code\tnote\n";
    CASES.forEach(c=>{const d=dec[c.token]||{};
      s+=[c.token,d.decision||"",d.reason||"",(d.note||"").replace(/[\t\n\r]/g," ")].join("\t")+"\n";});
    blob=new Blob([s],{type:"text/tab-separated-values"});name="NMV1C_ADJUDICATED_120.tsv";}
  else{const o=CASES.map(c=>({token:c.token,decision:(dec[c.token]||{}).decision||"",
      reason_code:(dec[c.token]||{}).reason||"",note:(dec[c.token]||{}).note||""}));
    blob=new Blob([JSON.stringify({exported:new Date().toISOString(),decisions:o},null,1)],
      {type:"application/json"});name="NMV1C_ADJUDICATED_120.json";}
  const u=URL.createObjectURL(blob),a=document.createElement("a");
  a.href=u;a.download=name;document.body.appendChild(a);a.click();
  setTimeout(()=>{URL.revokeObjectURL(u);a.remove();},0);}
function banner(){
  const b=document.getElementById("nostore");if(!b)return;
  b.style.display=STORE_OK?"none":"block";
  if(!STORE_OK)b.innerHTML="<b>Autosave is unavailable in this browser</b> ("+esc(STORE_WHY)+
    "). Your decisions are held in memory only and will be lost if this page is reloaded or "+
    "closed. Export the JSON regularly, and use <b>Restore progress</b> to load it back. "+
    "Opening this file directly from disk in Chrome, Edge or Firefox normally restores autosave.";}
function imp(f){const r=new FileReader();
  r.onload=()=>{let o;try{o=JSON.parse(r.result);}catch(e){alert("That file is not readable JSON.");return;}
    const rows=(o&&o.decisions)||[];let n=0;
    rows.forEach(d=>{if(d&&d.token){dec[d.token]={decision:d.decision||"",
      reason:d.reason_code||"",note:d.note||""};n++;}});
    save();render();alert("Restored "+n+" case(s).");};
  r.readAsText(f);}
try{render();}
catch(e){document.getElementById("case").innerHTML=
  "<p class='warn'>The application could not start: "+esc(e&&e.message||e)+"</p>";}
banner();
""")
    H.append("</script></html>")
    open(app, "w", encoding="utf-8", newline="\n").write("\n".join(H))

    # ---------------- manifest and sealed key ----------------
    man = os.path.join(D, "NMV1C_SAMPLE_MANIFEST%s.json" % SFX)
    MM = {"manifest": "NMV1C_SAMPLE_MANIFEST", "builder": VERSION,
          "frozen_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
          "protocol_sha256": PROTO_SHA, "seed": SEED, "n": len(fresh),
          "strata": {st: {"pool": len([b for b in avail if gt[b]["rule_id"] == st]),
                          "allocated": ALLOC[st], "drawn": len(sel[st]),
                          "weight": len([b for b in avail if gt[b]["rule_id"] == st]) / ALLOC[st],
                          "token_hash": hashlib.sha256("|".join(sorted(sel[st])).encode()).hexdigest()}
                     for st in ("A", "B", "D", "E")},
          "all_token_hash": hashlib.sha256("|".join(fresh).encode()).hexdigest(),
          "overlap_with_previous_packages": 0,
          "states_sampled": ["MOBILE", "QUIESCENT"],
          "states_not_sampled": ["NON_EVALUABLE"],
          "method_x_removed": True}
    json.dump(MM, open(man, "w", encoding="utf-8", newline="\n"), indent=2)
    key = os.path.join(D, "NMV1C_UNBLINDING_KEY%s.tsv" % SFX)
    with open(key, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("token\tblock_id\tstratum\trule_id\tmachine_state\tweight\n")
        for b in dsort(fresh, SEED, "nmv1ctoken"):
            rid = gt[b]["rule_id"]
            ms = "MOBILE" if rid in "ABCD" else ("QUIESCENT" if rid == "E" else "NON_EVALUABLE")
            fh.write("%s\t%s\t%s\t%s\t%s\t%.6f\n" % (tok[b], b, rid, rid, ms, weight[b]))
    print("\n  app      : %s" % sha(app))
    print("  manifest : %s" % sha(man))
    print("  key      : %s  (sealed)" % sha(key))


if __name__ == "__main__":
    main()
