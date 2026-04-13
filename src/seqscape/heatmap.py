from __future__ import annotations

import html
import sys
from pathlib import Path
from typing import List, Optional, Sequence


def maybe_heatmap(
    matrix: List[List[Optional[float]]],
    row_labels: List[str],
    col_labels: List[str],
    out_path: Path,
    title: str,
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"[warn] heatmap skipped ({out_path.name}): {exc}", file=sys.stderr)
        return

    n_rows = max(1, len(row_labels))
    n_cols = max(1, len(col_labels))
    fig_w = min(40.0, max(8.0, 0.28 * n_cols))
    fig_h = min(30.0, max(6.0, 0.34 * n_rows))

    m = [[float("nan") if v is None else float(v) for v in row] for row in matrix]
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad(color="#d9d9d9")
    cmap.set_under(color="#d9d9d9")
    im = ax.imshow(m, aspect="auto", interpolation="nearest", vmin=70, vmax=100, cmap=cmap)
    ax.set_title(title)
    ax.set_ylabel("Identified")
    ax.set_xlabel("Reference")
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=7)
    ax.set_xticks(range(len(col_labels)))
    if len(col_labels) <= 80:
        ax.set_xticklabels(col_labels, rotation=90, fontsize=6)
    else:
        ax.set_xticklabels([])
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("% similarity")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)


def color_for_value(v: Optional[float]) -> str:
    if v is None or v < 70.0:
        return "#d9d9d9"
    x = (max(70.0, min(100.0, v)) - 70.0) / 30.0
    r0, g0, b0 = (247, 251, 255)
    r1, g1, b1 = (8, 81, 156)
    r = int(r0 + (r1 - r0) * x)
    g = int(g0 + (g1 - g0) * x)
    b = int(b0 + (b1 - b0) * x)
    return f"#{r:02x}{g:02x}{b:02x}"


def write_html_heatmap(
    path: Path,
    row_labels: Sequence[str],
    col_labels: Sequence[str],
    matrix: List[List[Optional[float]]],
    title: str,
    filter_matrix: Optional[List[List[Optional[float]]]] = None,
    filter_label: str = "",
    enable_row_select: bool = False,
    col_labels_short: Optional[Sequence[str]] = None,
    row_export_ids: Optional[Sequence[str]] = None,
    export_filename: str = "selected_rows.txt",
    row_axis_label: str = "identified genomes",
    col_axis_label: str = "references",
) -> None:
    def esc(x: object) -> str:
        return html.escape(str(x), quote=True)

    has_short_labels = bool(col_labels_short) and (len(col_labels_short) == len(col_labels))
    has_row_export_ids = bool(row_export_ids) and (len(row_export_ids) == len(row_labels))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write("<!doctype html>\n")
        fh.write("<html><head><meta charset='utf-8'>\n")
        fh.write("<title>{}</title>\n".format(esc(title)))
        fh.write("<style>\n")
        fh.write("body { font-family: -apple-system, BlinkMacSystemFont, Segoe UI, sans-serif; margin: 12px; }\n")
        fh.write("h1 { font-size: 16px; margin: 0 0 10px 0; }\n")
        fh.write(".meta { margin: 0 0 10px 0; color: #444; font-size: 12px; }\n")
        fh.write(".controls { display: flex; flex-wrap: wrap; gap: 10px 16px; align-items: center; margin: 0 0 10px 0; font-size: 12px; }\n")
        fh.write(".controls label { display: inline-flex; gap: 6px; align-items: center; }\n")
        fh.write(".controls input[type='number'] { width: 72px; }\n")
        fh.write(".controls button { font-size: 12px; padding: 2px 8px; }\n")
        fh.write(".wrap { border: 1px solid #cfd8dc; max-height: 84vh; overflow: auto; }\n")
        fh.write("table { border-collapse: collapse; width: max-content; }\n")
        fh.write("th, td { border: 1px solid #e0e0e0; padding: 2px 6px; font-size: 11px; white-space: nowrap; text-align: center; }\n")
        fh.write("thead th { position: sticky; top: 0; background: #f8fafc; z-index: 3; }\n")
        fh.write("thead th.colhdr { height: 190px; min-width: 28px; padding: 4px 2px; vertical-align: bottom; }\n")
        fh.write("thead th.colhdr > span { display: inline-block; writing-mode: vertical-rl; transform: rotate(180deg); }\n")
        fh.write("body.hide-ref-labels thead th.colhdr > span { visibility: hidden; }\n")
        fh.write("tbody th { position: sticky; left: 0; background: #ffffff; z-index: 2; text-align: left; }\n")
        fh.write("tbody th.rowhdr label { display: inline-flex; gap: 6px; align-items: center; }\n")
        fh.write("thead th.corner { left: 0; z-index: 4; }\n")
        fh.write("</style></head><body>\n")
        fh.write(f"<h1>{esc(title)}</h1>\n")
        fh.write(f"<div class='meta'>Rows: {esc(row_axis_label)} | Columns: {esc(col_axis_label)} | Scroll down to keep column labels fixed.</div>\n")
        if filter_matrix is not None or enable_row_select or has_short_labels:
            fh.write("<div class='controls'>\n")
            if has_short_labels:
                fh.write("<button id='toggleNameModeBtn' type='button' data-mode='long'>Use Short Names</button>\n")
                fh.write("<label><input id='showRefLabels' type='checkbox' checked> show references</label>\n")
            if filter_matrix is not None:
                flab = esc(filter_label if filter_label else "qcov")
                fh.write(f"<label>{flab} min: <input id='qcovMinRange' type='range' min='0' max='100' step='0.1' value='0'><input id='qcovMinNumber' type='number' min='0' max='100' step='0.1' value='0'></label>\n")
                fh.write("<label><input id='hideNoPassRows' type='checkbox'> hide rows with no passing cells</label>\n")
                fh.write("<label><input id='hideNoPassCols' type='checkbox'> hide columns with no passing cells</label>\n")
            if enable_row_select:
                fh.write("<button id='selectVisibleBtn' type='button'>Select Visible</button>\n")
                fh.write("<button id='selectAllBtn' type='button'>Select All</button>\n")
                fh.write("<button id='clearSelectionBtn' type='button'>Clear</button>\n")
                fh.write("<button id='exportSelectionBtn' type='button'>Export</button>\n")
                fh.write("<span>selected: <strong id='selectedCount'>0</strong></span>\n")
            fh.write("</div>\n")
        fh.write("<div class='wrap'><table>\n<thead><tr><th class='corner'>identified_id</th>")
        for idx, c in enumerate(col_labels):
            c_short = col_labels_short[idx] if has_short_labels and col_labels_short else c
            fh.write(f"<th class='colhdr'><span class='coltxt' data-long='{esc(c)}' data-short='{esc(c_short)}'>{esc(c)}</span></th>")
        fh.write("</tr></thead>\n<tbody>\n")
        for ridx, rlab in enumerate(row_labels):
            fh.write("<tr>")
            erow = esc(rlab)
            export_id = row_export_ids[ridx] if has_row_export_ids and row_export_ids is not None else rlab
            eexport = esc(export_id)
            if enable_row_select:
                fh.write(f"<th class='rowhdr'><label><input class='row-select' type='checkbox' data-row='{eexport}'><span>{erow}</span></label></th>")
            else:
                fh.write(f"<th>{erow}</th>")
            row = matrix[ridx] if ridx < len(matrix) else []
            frow = filter_matrix[ridx] if (filter_matrix is not None and ridx < len(filter_matrix)) else []
            for cidx, _ in enumerate(col_labels):
                v = row[cidx] if cidx < len(row) else None
                fv = frow[cidx] if cidx < len(frow) else None
                color = color_for_value(v)
                txt = "NA" if v is None else f"{v:.1f}"
                qtxt = "" if fv is None else f"{fv:.3f}"
                fh.write(f"<td class='hmcell' data-missing='{1 if v is None else 0}' data-base-color='{esc(color)}' data-base-text='{esc(txt)}' data-qcov='{esc(qtxt)}' style='background:{esc(color)}'>{esc(txt)}</td>")
            fh.write("</tr>\n")
        fh.write("</tbody></table></div>\n")
        fh.write("<script>\n(function(){\n")
        fh.write("const qRange=document.getElementById('qcovMinRange'); const qNumber=document.getElementById('qcovMinNumber'); const hideNoPassRows=document.getElementById('hideNoPassRows'); const hideNoPassCols=document.getElementById('hideNoPassCols'); const selectVisibleBtn=document.getElementById('selectVisibleBtn'); const selectAllBtn=document.getElementById('selectAllBtn'); const clearSelectionBtn=document.getElementById('clearSelectionBtn'); const exportSelectionBtn=document.getElementById('exportSelectionBtn'); const toggleNameModeBtn=document.getElementById('toggleNameModeBtn'); const showRefLabels=document.getElementById('showRefLabels'); const selectedCount=document.getElementById('selectedCount'); const rowChecks=Array.from(document.querySelectorAll('.row-select')); const colHeaders=Array.from(document.querySelectorAll('thead th.colhdr')); const colTexts=Array.from(document.querySelectorAll('span.coltxt'));\n")
        fh.write("function applyNameMode(){ if(!toggleNameModeBtn) return; const mode=toggleNameModeBtn.dataset.mode||'long'; for(const sp of colTexts){ const longName=sp.dataset.long||''; const shortName=sp.dataset.short||longName; sp.textContent=(mode==='short')?shortName:longName; }}\n")
        fh.write("function isVisibleRowCheckbox(ch){ const tr=ch.closest('tr'); if(!tr) return false; return tr.style.display !== 'none'; }\n")
        fh.write("function updateSelectedCount(){ if(!selectedCount) return; let n=0; for(const ch of rowChecks) if(ch.checked) n += 1; selectedCount.textContent=String(n); }\n")
        fh.write("function applyQcovFilter(){ const minQ=qRange?parseFloat(qRange.value||'0'):0; const hideRows=hideNoPassRows&&hideNoPassRows.checked&&!!qRange; const hideCols=hideNoPassCols&&hideNoPassCols.checked&&!!qRange; const rows=Array.from(document.querySelectorAll('tbody tr')); const colPass=new Array(colHeaders.length).fill(false); for(const tr of rows){ let rowPass=false; const cells=Array.from(tr.querySelectorAll('td.hmcell')); for(let ci=0; ci<cells.length; ci+=1){ const td=cells[ci]; const miss=td.dataset.missing==='1'; const qraw=td.dataset.qcov||''; const q=qraw===''?NaN:parseFloat(qraw); let show=true; if(miss){ show=false; } else if(qRange){ show=!Number.isNaN(q) && q>=minQ; } if(show){ td.style.background=td.dataset.baseColor||'#ffffff'; td.textContent=td.dataset.baseText||''; rowPass=true; if(ci<colPass.length) colPass[ci]=true; } else { td.style.background='#d9d9d9'; td.textContent=td.dataset.baseText||'NA'; }} tr.style.display=hideRows && !rowPass ? 'none' : ''; } for(let ci=0; ci<colHeaders.length; ci+=1){ const showCol=!hideCols || colPass[ci]; colHeaders[ci].style.display=showCol ? '' : 'none'; } for(const tr of rows){ const cells=Array.from(tr.querySelectorAll('td.hmcell')); for(let ci=0; ci<cells.length; ci+=1){ const showCol=!hideCols || colPass[ci]; cells[ci].style.display=showCol ? '' : 'none'; }} updateSelectedCount(); }\n")
        fh.write("if(qRange&&qNumber){ qRange.addEventListener('input', function(){ qNumber.value=qRange.value; applyQcovFilter(); }); qNumber.addEventListener('input', function(){ let v=parseFloat(qNumber.value||'0'); if(Number.isNaN(v)) v=0; if(v<0) v=0; if(v>100) v=100; qRange.value=String(v); qNumber.value=String(v); applyQcovFilter(); }); }\n")
        fh.write("if(hideNoPassRows) hideNoPassRows.addEventListener('change', applyQcovFilter); if(hideNoPassCols) hideNoPassCols.addEventListener('change', applyQcovFilter);\n")
        fh.write("if(toggleNameModeBtn) toggleNameModeBtn.addEventListener('click', function(){ const cur=toggleNameModeBtn.dataset.mode||'long'; if(cur==='long'){ toggleNameModeBtn.dataset.mode='short'; toggleNameModeBtn.textContent='Use Long Names'; } else { toggleNameModeBtn.dataset.mode='long'; toggleNameModeBtn.textContent='Use Short Names'; } applyNameMode(); });\n")
        fh.write("if(showRefLabels) showRefLabels.addEventListener('change', function(){ document.body.classList.toggle('hide-ref-labels', !showRefLabels.checked); });\n")
        fh.write("if(selectVisibleBtn) selectVisibleBtn.addEventListener('click', function(){ for(const ch of rowChecks){ if(isVisibleRowCheckbox(ch)) ch.checked=true; } updateSelectedCount(); });\n")
        fh.write("if(selectAllBtn) selectAllBtn.addEventListener('click', function(){ for(const ch of rowChecks) ch.checked=true; updateSelectedCount(); });\n")
        fh.write("if(clearSelectionBtn) clearSelectionBtn.addEventListener('click', function(){ for(const ch of rowChecks) ch.checked=false; updateSelectedCount(); });\n")
        fh.write("for(const ch of rowChecks) ch.addEventListener('change', updateSelectedCount);\n")
        fh.write(f"if(exportSelectionBtn) exportSelectionBtn.addEventListener('click', function(){{ const names=rowChecks.filter(ch => ch.checked).map(ch => ch.dataset.row || ''); if(!names.length){{ alert('No rows selected.'); return; }} const blob=new Blob([names.join('\\\\n') + '\\\\n'], {{type:'text/plain;charset=utf-8'}}); const a=document.createElement('a'); const url=URL.createObjectURL(blob); a.href=url; a.download={export_filename!r}; document.body.appendChild(a); a.click(); document.body.removeChild(a); setTimeout(function(){{ URL.revokeObjectURL(url); }}, 1000); }});\n")
        fh.write("applyNameMode(); applyQcovFilter();})();\n</script>\n</body></html>\n")
