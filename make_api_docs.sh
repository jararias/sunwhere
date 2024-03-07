# python3 -m pip install https://github.com/jararias/lazydocs

# lazydocs sunwhere --overview-file index.md --output-path api-docs/
lazydocs sunwhere/base.py sunwhere/usecases.py --src-base-url https://github.com/jararias/sunwhere/blob/master/ --output-path api-docs/
# rm -f api-docs/.pages api-docs/core.md api-docs/libspa.md
# rm -r api-docs/utils.md api-docs/version.md
mv api-docs/base.py.md api-docs/base.md
mv api-docs/usecases.py.md api-docs/usecases.md