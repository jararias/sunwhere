# cd <path to sunwhere>
coverage run -m pytest -v --junit-xml=junit.xml tests/
coverage xml
genbadge tests -i junit.xml -o assets/tests-badge.svg
genbadge coverage -i coverage.xml -o assets/coverage-badge.svg
rm -f junit.xml coverage.xml
