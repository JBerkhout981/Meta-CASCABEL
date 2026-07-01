import os
from datetime import datetime

def read_file(path, max_lines=200):
    if not os.path.exists(path):
        return f"<p style='color:red;'>Missing: {path}</p>"

    with open(path) as f:
        lines = f.readlines()[:max_lines]
    return "<pre>" + "".join(lines) + "</pre>"


out = snakemake.output.html

sections = []

sections.append(f"<h1>Metagenomics report</h1>")
sections.append(f"<p>Generated: {datetime.now()}</p>")

# Mapping
sections.append("<h2>Mapping / Flagstat</h2>")
sections.append(read_file(snakemake.input.flagstat))

# Taxonomy
sections.append("<h2>Taxonomy</h2>")
sections.append(read_file(snakemake.input.taxonomy))

# Binning
sections.append("<h2>Binning summary</h2>")
sections.append(read_file(snakemake.input.binning_summary))

sections.append("<h2>Abundance</h2>")
sections.append(read_file(snakemake.input.abundance))

# Unbinned
sections.append("<h2>Unbinned</h2>")
sections.append(read_file(snakemake.input.unbinned))

# Final bins
sections.append("<h2>Final bins</h2>")
sections.append(read_file(snakemake.input.final_bins))

# Coverage
sections.append("<h2>Coverage</h2>")
sections.append(read_file(snakemake.input.coverage))

html = """
<html>
<head>
<style>
body { font-family: Arial; margin: 40px; }
pre { background: #f4f4f4; padding: 10px; overflow-x: auto; }
h2 { border-bottom: 1px solid #ddd; padding-bottom: 5px; }
</style>
</head>
<body>
""" + "\n".join(sections) + "</body></html>"

with open(out, "w") as f:
    f.write(html)
