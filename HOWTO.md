
## Checklist to publish

1. Check that the **documentation** is published and works properly
2. Check all the **metadata** in `pyproject.toml`, in particular, that the Documentation entry under `[project.urls]` in `pyproject.toml` points out to the documentation site.
3. Check that the **version number** is correct. For a version number *x.y.z*, you must upgrade *z* for minor updates, *y* for more important changes that maintain backward compatibility and *x* for significant changes that can affect the backward compatibility. Whenever you upgrade *y* or *x* the numbers to the right (that is, *z* and *y.z*, respectively) must be set to *0*.
4. **Build** the package:
    ```
    uv build  # it must create dist/sunwhere-x.y.z* (.whl & .tar.gz)
    ```
5. **Publish** the package **in test.pypi**:
    ```
    uv publish --index-url https://test.pypi.org/simple --token <test-pypi-token>
    ```
   And then check that it installs properly:
    ```bash
    uv init --python 3.xy dummy-env
    cd dummy-env
    uv add --index-url https://test.pypi.org/simple/ sunwhere
    # check that it works as expected...
    # ...
    cd .. && rm -rf dummy-env
    ```
7. If all goes well, **publish in the primary repository** (pypi.org):
    ```bash
    uv publish --token <pypi-token>
    # check that it installs properly...
    uv init --python 3.xy dummy-env
    cd dummy-env
    uv add sunwhere
    # and that it works as expected...
    # ...
    cd .. && rm -rf dummy-env
    ```
8. Create a **GitHub Release**:
    ```bash
    git tag -a vx.y.z -m "Release vx.y.z"
    git push origin vx.y.z
    ```
    Then go to *GitHub* -> *Releases* -> *Draft a new release* and link the tag `vx.y.z` to the release.
