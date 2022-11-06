## Deploy

```bash
# Documentation
# pip install sphinx nbsphinx myst_parser
sphinx-build -b html . build/html

# Pypi package
rm -r dist
# change version in setup.py
# change something in code
python3 -m build
python3 -m twine upload --repository pypi dist/*
Uploading distributions to https://upload.pypi.org/legacy/
# Enter your username: ...
# Enter your password: ...
```
