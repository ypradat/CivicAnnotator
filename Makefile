# Makefile macros (or variables)
PYTHON ?= python
PYTEST ?= pytest
CTAGS ?= ctags

# .PHONY defines parts of the makefile that are not dependant on any specific file
# This is most often used to store functions
.PHONY = clean develop install tags test

clean:
	rm -f tags
	rm -rf build
	rm -rf dist
	rm -rf *.egg-info

tags:
	$(CTAGS) --python-kinds=-i --exclude=*/tests/* -R .

test:
	@echo "---------------Run prettypy tests-----------------"
	$(PYTHON) -m pip install pytest pytest-cov
	$(PYTEST) --cov-config=.coveragerc --cov-report xml:coverage.xml --cov . tests
