PY ?= python3

.PHONY: install test test-all artifacts verify zeros clean dist

install:
	$(PY) -m pip install -e ".[test]"

test:
	$(PY) -m pytest -q -m "not slow"

test-all:
	$(PY) -m pytest -q

artifacts:
	$(PY) scripts/generate_artifacts.py --out artifacts

verify:
	$(PY) -m hausdorff_certificates.verify artifacts/*.cert.json

zeros:
	rm -f data/zeta_zeros_60_dps60.json
	$(PY) -c "from hausdorff_certificates.zeta import load_or_compute_zeros as f; f('data/zeta_zeros_60_dps60.json', 60, dps=60); print('zeros regenerated')"

clean:
	rm -rf artifacts_ci build dist *.egg-info src/*.egg-info .pytest_cache
	find . -name __pycache__ -type d -exec rm -rf {} +

dist: clean
	git archive --format=zip -o hausdorff-certificates.zip HEAD
