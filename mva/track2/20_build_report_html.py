"""Render the Track 2 report markdown as a standalone HTML page.

Content fidelity matters more than cleverness here: the markdown is the
submission, so this converts every block rather than hand-copying prose that
could silently drop a caveat. Design lives in the template below; this half
only turns markdown constructs into semantic HTML.

Run from the repo root: python3 -I mva/track2/20_build_report_html.py
"""
import html
import os
import re
import sys

SRC = "reports/track2-mechanism-and-repurposing.md"
OUT = "reports/track2-report.html"

if not os.path.exists(SRC):
    sys.exit(f"{SRC} not found - run from the repository root")

# ---------------------------------------------------------------- inline
def inline(text):
    """Markdown inline -> HTML, with code spans and links protected from the
    emphasis pass so asterisks inside them survive."""
    slots = []

    def stash(h):
        slots.append(h)
        return f"\x00{len(slots)-1}\x00"

    # code spans first: nothing inside them is markdown
    text = re.sub(r"`([^`]+)`",
                  lambda m: stash(f"<code>{html.escape(m.group(1))}</code>"),
                  text)
    # autolinks <https://...>
    text = re.sub(r"<(https?://[^>\s]+)>",
                  lambda m: stash(f'<a href="{html.escape(m.group(1))}">'
                                  f'{html.escape(m.group(1))}</a>'),
                  text)
    # [label](url) - label may itself contain emphasis, so recurse
    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)",
                  lambda m: stash(f'<a href="{html.escape(m.group(2))}">'
                                  f'{inline(m.group(1))}</a>'),
                  text)
    # escaped asterisk, e.g. UGT1A1\*28
    text = text.replace(r"\*", "\x01")

    text = html.escape(text, quote=False)
    text = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", text)
    text = re.sub(r"\*([^*]+)\*", r"<em>\1</em>", text)
    text = text.replace("\x01", "*")

    for i, h in enumerate(slots):
        text = text.replace(f"\x00{i}\x00", h)
    return text


def slug(heading):
    m = re.match(r"(\d+(?:\.\d+)?)\.?\s", heading)
    if m:
        return "s" + m.group(1).replace(".", "-")
    m = re.match(r"Appendix ([A-C])", heading)
    if m:
        return "app-" + m.group(1).lower()
    return re.sub(r"[^a-z0-9]+", "-", heading.lower()).strip("-")[:40]


def table_class(header_cells):
    first = re.sub(r"[^a-z#]", "", header_cells[0].lower())
    joined = " ".join(c.lower() for c in header_cells)
    if first == "#" and "candidate" in joined:
        return " scored"
    if "gate" in joined or "observed" in joined:
        return " decision"
    return ""


# ---------------------------------------------------------------- blocks
def convert(md):
    lines = md.split("\n")
    out, i = [], 0
    sections = []          # (id, label) for the nav
    open_section = False

    while i < len(lines):
        ln = lines[i]

        # fenced code
        if ln.startswith("```"):
            i += 1
            buf = []
            while i < len(lines) and not lines[i].startswith("```"):
                buf.append(lines[i]); i += 1
            i += 1
            out.append('<div class="scroll"><pre class="ascii">'
                       + html.escape("\n".join(buf)) + "</pre></div>")
            continue

        # headings
        m = re.match(r"^(#{1,3})\s+(.*)$", ln)
        if m:
            level, txt = len(m.group(1)), m.group(2).strip()
            if level == 1:
                i += 1
                continue                      # the hero carries the H1
            sid = slug(txt)
            if level == 2:
                if open_section:
                    out.append("</section>")
                out.append(f'<section id="{sid}">')
                open_section = True
                num = re.match(r"(\d+|Appendix [A-C])", txt)
                sections.append((sid, txt))
                label = re.sub(r"^(\d+)\.\s*", "", txt)
                out.append(f'<h2><span class="hn">{html.escape(num.group(1)) if num else ""}</span>'
                           f'{inline(label)}</h2>')
            else:
                label = re.sub(r"^(\d+\.\d+)\s*", "", txt)
                sub = re.match(r"(\d+\.\d+)", txt)
                out.append(f'<h3 id="{sid}">'
                           f'<span class="hn">{html.escape(sub.group(1)) if sub else ""}</span>'
                           f'{inline(label)}</h3>')
            i += 1
            continue

        # figure
        m = re.match(r"^!\[([^\]]*)\]\(([^)]+)\)\s*$", ln)
        if m:
            alt, src = html.escape(m.group(1)), html.escape(m.group(2))
            out.append(f'<figure><img src="{src}" alt="{alt}" loading="lazy">'
                       f'<figcaption>{alt}</figcaption></figure>')
            i += 1
            continue

        # horizontal rule
        if re.match(r"^---+\s*$", ln):
            i += 1
            continue

        # blockquote (may contain its own heading)
        if ln.startswith(">"):
            buf = []
            while i < len(lines) and (lines[i].startswith(">") or
                                      (buf and lines[i].strip() == ""
                                       and i + 1 < len(lines)
                                       and lines[i + 1].startswith(">"))):
                buf.append(re.sub(r"^>\s?", "", lines[i])); i += 1
            inner, cls = [], "note"
            for b in buf:
                if not b.strip():
                    continue
                hm = re.match(r"^#{1,4}\s+(.*)$", b)
                if hm:
                    cls = "note warn"
                    inner.append(f"<p class='note-h'>{inline(hm.group(1))}</p>")
                else:
                    inner.append(f"<p>{inline(b)}</p>")
            out.append(f'<aside class="{cls}">' + "".join(inner) + "</aside>")
            continue

        # table
        if ln.startswith("|") and i + 1 < len(lines) and \
                re.match(r"^\|[\s:|-]+\|\s*$", lines[i + 1]):
            def cells(row):
                return [c.strip() for c in row.strip().strip("|").split("|")]
            head = cells(ln)
            i += 2
            body = []
            while i < len(lines) and lines[i].startswith("|"):
                body.append(cells(lines[i])); i += 1
            cls = table_class(head)
            t = [f'<div class="scroll"><table class="tbl{cls}"><thead><tr>']
            for c in head:
                t.append(f"<th>{inline(c)}</th>")
            t.append("</tr></thead><tbody>")
            for row in body:
                t.append("<tr>" + "".join(f"<td>{inline(c)}</td>" for c in row)
                         + "</tr>")
            t.append("</tbody></table></div>")
            out.append("".join(t))
            continue

        # ordered list
        if re.match(r"^\d+\.\s+", ln):
            items = []
            while i < len(lines) and re.match(r"^\d+\.\s+", lines[i]):
                items.append(re.sub(r"^\d+\.\s+", "", lines[i])); i += 1
            out.append("<ol>" + "".join(f"<li>{inline(x)}</li>" for x in items)
                       + "</ol>")
            continue

        # unordered list
        if re.match(r"^[-*]\s+", ln):
            items = []
            while i < len(lines) and re.match(r"^[-*]\s+", lines[i]):
                items.append(re.sub(r"^[-*]\s+", "", lines[i])); i += 1
            cls = ' class="refs"' if len(items) > 12 else ""
            out.append(f"<ul{cls}>"
                       + "".join(f"<li>{inline(x)}</li>" for x in items)
                       + "</ul>")
            continue

        # paragraph. NOTE: a "|" or "![" line that failed its block test
        # above would otherwise consume nothing and spin here forever, so the
        # loop below is guarded to always advance at least one line.
        if ln.strip():
            buf = [ln]
            i += 1
            while i < len(lines) and lines[i].strip() and \
                    not lines[i].startswith(("|", ">", "#", "```", "![")) and \
                    not re.match(r"^(\d+\.|[-*])\s+", lines[i]):
                buf.append(lines[i]); i += 1
            out.append("<p>" + inline(" ".join(buf)) + "</p>")
            continue

        i += 1

    if open_section:
        out.append("</section>")
    return "\n".join(out), sections


md = open(SRC).read()
# the hero is hand-built from the front matter; drop it from the body
body_md = md.split("---\n", 1)[1] if "---\n" in md else md
body_html, sections = convert(body_md)

nav = "".join(
    f'<a href="#{sid}">{html.escape(re.sub(r"^(\d+)\.\s*", r"\1 · ", t))}</a>'
    for sid, t in sections)

TEMPLATE = """<title>BubR1 Repurposing Dossier</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:ital,wght@0,400;0,500;0,600;1,400&family=IBM+Plex+Mono:wght@400;500&family=Source+Serif+4:ital,opsz,wght@0,8..60,400;0,8..60,600;1,8..60,400&display=swap">
<style>
/* Layout: one reading column with a sticky section rail; figures and tables
   break wider than the measure. Palette keyed to the committed figure PNGs. */
:root {
  --paper:#fbfbfa; --panel:#f3f3f1; --rule:#e0e0dc;
  --ink:#1a1b1a; --ink-2:#55574f; --ink-3:#83857c;
  --accent:#1f5fab; --accent-soft:#e8eff8;
  --warn:#c2521f; --warn-soft:#fbeadf;
  --ok:#17795a;
  --sans:"IBM Plex Sans",system-ui,-apple-system,"Segoe UI",sans-serif;
  --serif:"Source Serif 4",Georgia,"Times New Roman",serif;
  --mono:"IBM Plex Mono",ui-monospace,"SF Mono",Menlo,monospace;
  --measure:70ch;
  color-scheme:light;
}
@media (prefers-color-scheme:dark){ :root:not([data-theme="light"]){
  --paper:#15161a; --panel:#1d1f24; --rule:#30333a;
  --ink:#e9eaec; --ink-2:#a9adb6; --ink-3:#7d818b;
  --accent:#6fa8ea; --accent-soft:#1a2738;
  --warn:#f08a56; --warn-soft:#2e2019;
  --ok:#49c49b;
  color-scheme:dark;
}}
:root[data-theme="dark"]{
  --paper:#15161a; --panel:#1d1f24; --rule:#30333a;
  --ink:#e9eaec; --ink-2:#a9adb6; --ink-3:#7d818b;
  --accent:#6fa8ea; --accent-soft:#1a2738;
  --warn:#f08a56; --warn-soft:#2e2019;
  --ok:#49c49b;
  color-scheme:dark;
}
*{box-sizing:border-box}
body{background:var(--paper);color:var(--ink);font-family:var(--serif);
  font-size:17px;line-height:1.62;-webkit-text-size-adjust:100%}
body::before{content:"";position:fixed;inset:0 0 auto;
  height:env(safe-area-inset-top,0px);background:var(--paper);z-index:20}
.gutter{padding-left:max(18px,env(safe-area-inset-left,0px));
  padding-right:max(18px,env(safe-area-inset-right,0px))}
.wrap{max-width:var(--measure);margin:0 auto}
a{color:var(--accent);text-underline-offset:2px;
  text-decoration-color:color-mix(in srgb,var(--accent) 35%,transparent)}
a:hover{text-decoration-color:var(--accent)}
:focus-visible{outline:2px solid var(--accent);outline-offset:3px;border-radius:2px}
code{font-family:var(--mono);font-size:.86em;background:var(--panel);
  padding:.08em .34em;border-radius:3px;color:var(--ink-2);
  word-break:break-word}
strong{font-weight:600}

/* ---- masthead ---- */
.hero{border-bottom:1px solid var(--rule);padding-block:34px 26px}
.eyebrow{font-family:var(--sans);font-size:11.5px;font-weight:600;
  letter-spacing:.1em;text-transform:uppercase;color:var(--ink-3);
  margin:0 0 14px}
h1{font-family:var(--sans);font-weight:600;font-size:clamp(28px,6.2vw,42px);
  line-height:1.1;letter-spacing:-.02em;margin:0;text-wrap:balance}
.dek{font-size:clamp(17px,2.6vw,19.5px);color:var(--ink-2);margin:14px 0 0;
  max-width:62ch}
.meta{display:flex;flex-wrap:wrap;gap:6px 16px;margin-top:20px;
  font-family:var(--sans);font-size:12.5px;color:var(--ink-3)}
.meta b{font-weight:500;color:var(--ink-2)}
.facts{display:grid;grid-template-columns:repeat(auto-fit,minmax(118px,1fr));
  gap:14px 20px;margin-top:26px;padding-top:20px;border-top:1px solid var(--rule)}
.fact .n{font-family:var(--sans);font-weight:600;font-size:21px;
  font-variant-numeric:tabular-nums;letter-spacing:-.01em;display:block}
.fact .l{font-family:var(--sans);font-size:11.5px;line-height:1.35;
  color:var(--ink-3);display:block;margin-top:3px}
.fact.neg .n{color:var(--warn)}

/* ---- section rail ---- */
nav{position:sticky;top:env(safe-area-inset-top,0px);z-index:19;
  background:color-mix(in srgb,var(--paper) 92%,transparent);
  backdrop-filter:saturate(1.6) blur(8px);
  border-bottom:1px solid var(--rule)}
.navscroll{display:flex;gap:2px;overflow-x:auto;scrollbar-width:none;
  padding-block:9px;max-width:var(--measure);margin:0 auto}
.navscroll::-webkit-scrollbar{display:none}
nav a{font-family:var(--sans);font-size:12px;font-weight:500;
  white-space:nowrap;padding:5px 9px;border-radius:5px;color:var(--ink-2);
  text-decoration:none}
nav a:hover{background:var(--panel);color:var(--ink)}

/* ---- prose ---- */
main{padding-block:8px 60px}
section{padding-block:34px 4px;border-top:1px solid var(--rule)}
section:first-of-type{border-top:0}
h2{font-family:var(--sans);font-weight:600;font-size:clamp(21px,3.4vw,25px);
  line-height:1.22;letter-spacing:-.015em;margin:0 0 18px;text-wrap:balance;
  scroll-margin-top:72px}
h3{font-family:var(--sans);font-weight:600;font-size:17.5px;line-height:1.3;
  margin:36px 0 14px;text-wrap:balance;scroll-margin-top:72px}
.hn{font-family:var(--mono);font-weight:400;color:var(--accent);
  font-size:.78em;margin-right:.7em;letter-spacing:0}
h2 .hn:empty,h3 .hn:empty{display:none}
p{margin:0 0 15px}
ol,ul{margin:0 0 16px;padding-left:1.35em}
li{margin-bottom:9px}
ol li::marker{font-family:var(--sans);font-weight:600;color:var(--accent);
  font-size:.9em}
ul.refs{list-style:none;padding-left:0;font-size:14.5px;
  column-gap:34px;color:var(--ink-2)}
ul.refs li{margin-bottom:11px;padding-left:1.1em;text-indent:-1.1em;
  line-height:1.45;overflow-wrap:break-word}

/* ---- callouts ---- */
aside.note{background:var(--panel);border-left:3px solid var(--accent);
  border-radius:0 6px 6px 0;padding:15px 18px;margin:0 0 18px}
aside.note p{margin:0 0 9px}
aside.note p:last-child{margin:0}
aside.note.warn{background:var(--warn-soft);border-left-color:var(--warn)}
.note-h{font-family:var(--sans);font-weight:600;font-size:15.5px}

/* ---- figures ---- */
figure{margin:26px 0 24px}
figure img{display:block;width:100%;max-width:100%;height:auto;
  border:1px solid var(--rule);border-radius:7px;background:#fcfcfb}
figcaption{font-family:var(--sans);font-size:12px;color:var(--ink-3);
  margin-top:9px;line-height:1.4}

/* ---- tables ---- */
.scroll{overflow-x:auto;-webkit-overflow-scrolling:touch;margin:0 0 20px;
  border:1px solid var(--rule);border-radius:7px;min-width:0}
table.tbl{border-collapse:collapse;width:100%;min-width:560px;
  font-family:var(--sans);font-size:13.2px;line-height:1.48}
table.tbl th{text-align:left;font-weight:600;font-size:11.5px;
  letter-spacing:.045em;text-transform:uppercase;color:var(--ink-3);
  background:var(--panel);padding:9px 12px;border-bottom:1px solid var(--rule);
  white-space:nowrap;vertical-align:bottom}
table.tbl td{padding:11px 12px;border-bottom:1px solid var(--rule);
  vertical-align:top;color:var(--ink-2)}
table.tbl tbody tr:last-child td{border-bottom:0}
table.tbl td strong{color:var(--ink)}
table.tbl td code{font-size:.84em}
table.scored td:nth-child(5),table.scored td:nth-child(6),
table.scored td:nth-child(7){text-align:center;font-variant-numeric:tabular-nums;
  font-weight:600;color:var(--ink);white-space:nowrap}
table.scored td:nth-child(1),table.scored td:nth-child(3){text-align:center;
  font-weight:600;color:var(--ink);white-space:nowrap}
table.decision td:first-child{font-weight:600;color:var(--ink);
  white-space:nowrap}
pre.ascii{font-family:var(--mono);font-size:12px;line-height:1.5;margin:0;
  padding:16px;color:var(--ink-2);white-space:pre;tab-size:2}

/* ---- footer ---- */
footer{border-top:1px solid var(--rule);padding-block:24px 40px;
  font-family:var(--sans);font-size:12.5px;color:var(--ink-3)}
footer p{margin:0 0 7px;max-width:68ch}

@media (max-width:640px){
  body{font-size:16px}
  .facts{grid-template-columns:repeat(2,1fr)}
}
@media (prefers-reduced-motion:reduce){
  *{animation:none!important;transition:none!important;scroll-behavior:auto!important}
}
html{scroll-behavior:smooth}
</style>

<header class="hero gutter">
  <div class="wrap">
    <p class="eyebrow">Rare Disease, Real Kid · MVA Hackathon 2026 · Track 2</p>
    <h1>Mechanism and Drug Repurposing in <em>BUB1B</em>-associated MVA</h1>
    <p class="dek">Genotype-specific repurposing for a child with biallelic <em>BUB1B</em>: record-level assay evidence to choose candidates, and functional tests to separate restored accurate division from raised protein abundance, suppressed proliferation, or selective loss of the worst cells.</p>
    <div class="meta">
      <span>proband <b>PROBAND01</b></span>
      <span>Daniel Hwang · building on team <b>texdata</b></span>
      <span><b>24 October 2026</b></span>
      <span>code: <code>mva/track2/</code></span>
    </div>
    <div class="facts">
      <div class="fact"><span class="n">10</span><span class="l">experiments, E0&ndash;E9</span></div>
      <div class="fact"><span class="n">11</span><span class="l">candidates scored on four axes</span></div>
      <div class="fact"><span class="n">1,018</span><span class="l">ClinGen genes surveyed</span></div>
      <div class="fact"><span class="n">42:1</span><span class="l">inhibition to activation in ChEMBL</span></div>
      <div class="fact neg"><span class="n">0</span><span class="l">decision changes that survived adjudication</span></div>
    </div>
  </div>
</header>

<nav class="gutter"><div class="navscroll">__NAV__</div></nav>

<main class="gutter"><div class="wrap">
__BODY__
</div></main>

<footer class="gutter"><div class="wrap">
  <p><strong>Hypotheses for follow-up, not evidence that a medicine works.</strong> Nothing here has been tested in anyone with MVA1, and nothing substitutes for this child&rsquo;s oncology team.</p>
  <p>Written with Claude Opus 5 via Claude Code. Figures are generated from committed JSON by <code>17_figures.py</code> and <code>22_structure_figure.py</code>. Numbers in &sect;4, &sect;5 and the appendices trace to committed script outputs; &sect;3&rsquo;s graph-model figures do not, and say so. Corrections and their reasons are in the repository&rsquo;s commit history.</p>
</div></footer>
"""

page = TEMPLATE.replace("__NAV__", nav).replace("__BODY__", body_html)
open(OUT, "w").write(page)
print(f"wrote {OUT}  ({len(page):,} bytes, {len(sections)} sections)")
for sid, t in sections:
    print(f"  #{sid:<10} {t}")
