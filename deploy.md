## Deploy

- https://docs.python-guide.org/writing/documentation/
- https://docs.readthedocs.io/en/stable/tutorial/
- https://samnicholls.net/2016/06/15/how-to-sphinx-readthedocs/
- https://stackoverflow.com/questions/23211695/modifying-content-width-of-the-sphinx-theme-read-the-docs


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

Once the project is registered with readthedocs.org the documentation rebuilds each time there is a push to the Github repo (login to readthedocs to see the build progress/ any errors).

