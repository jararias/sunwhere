
.PHONY: ipython3
ipython3:
	@uv run ipython --profile-dir=.ipython

.PHONY: init-venv
init-venv:
	@uv sync --reinstall
	@attr -s com.dropbox.ignored -V 1 .venv  # instruct dropbox to ignore the .venv folder

.PHONY: clean
clean:
	@rm -rf build
	@rm -rf dist
	@rm -rf *.egg-info
	@rm -rf sunwhere/_version.py
	@rm -rf .pytest_cache
	@rm -f MANIFEST
	@find . -name "__pycache__" -print0 | xargs -0 -I {} /bin/rm -rf "{}"

.PHONY: tests
tests:
	@uv run coverage run -m pytest -v --junit-xml=junit.xml tests/
	@uv run coverage xml
	@uv run genbadge tests -i junit.xml -o assets/tests-badge.svg
	@uv run genbadge coverage -i coverage.xml -o assets/coverage-badge.svg
	@rm -f junit.xml coverage.xml
