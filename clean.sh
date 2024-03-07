rm -rf build
rm -rf dist
rm -rf *.egg-info
rm -rf sunwhere/_version.py
rm -rf .pytest_cache
rm -f MANIFEST
find . -name "__pycache__" -print0 | xargs -0 -I {} /bin/rm -rf "{}"
