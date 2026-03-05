
# Development

See the [tskit developer documentation](https://tskit.dev/tskit/docs/stable/development.html)
for the general development workflow (git, prek, testing, documentation).

Install development dependencies with:

```
$ uv sync
```

Run the unit tests with:

```
$ uv run pytest
```

You may need to regenerate some cached test fixtures occasionally (particularly
if getting cryptic errors when running the test suite). To do this, run:

```
rm -fR tests/data/cache/
```

and rerun tests as above.

## Debug utilities

The tree sequence files output during primary inference have a lot
of debugging metadata, and there are some developer tools for inspecting
this in the ``sc2ts.debug`` package. In particular, the ``ArgInfo``
class has a lot of useful utilities designed to be used in a Jupyter
notebook. Note that ``matplotlib`` is required for these. Use it like:

```python
import sc2ts.debug as sd
import tskit

ts = tskit.load("path_to_daily_inference.ts")
ai = sd.ArgInfo(ts)
ai # view summary in notebook
```


