"""
Sphinx extension that makes the plotly figures produced by the notebooks work in the book.

myst-nb can only render the text/html output of plotly figures, but the default
renderer depends on the environment (e.g. "vscode" when building from a VS Code
terminal) and may emit only application/vnd.plotly.v1+json, so the figures are
dropped from the book. The kernels inherit PLOTLY_RENDERER from the build process.
Export PLOTLY_RENDERER to override the default.

The HTML emitted by the plotly renderers also loads MathJax 2 from a CDN. This
replaces the MathJax 3/4 configuration of the book, so that equations are no longer
rendered, and the init cell imports an URL without extension that the CDN refuses.
Both scripts are removed from the cell outputs.
"""
import os
import re

from docutils import nodes

MATHJAX2_SCRIPT = re.compile(
    r'<script[^>]*src="[^"]*/mathjax/2\.[^"]*MathJax\.js[^"]*"[^>]*>\s*</script>'
)
PLOTLY_MODULE_IMPORT = re.compile(
    r'<script type="module">\s*import "https://cdn\.plot\.ly/plotly-[^"]*\.min"\s*</script>'
)


def strip_plotly_mathjax(app, doctree, docname):
    for node in doctree.findall(nodes.raw):
        if "html" not in node.get("format", "").split():
            continue
        text = node.astext()
        new = PLOTLY_MODULE_IMPORT.sub("", MATHJAX2_SCRIPT.sub("", text))
        if new != text:
            node.replace_self(nodes.raw(new, new, **node.attributes))


def setup(app):
    os.environ.setdefault("PLOTLY_RENDERER", "notebook_connected")
    app.connect("doctree-resolved", strip_plotly_mathjax)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
