

default: html_funcs

asciidoc_funcs:
	./extract_matlab_help.py ../msat/*.m

html_funcs: asciidoc_funcs
	(cd funcs ; for f in *.txt ; do asciidoc $$f ; done)

clean:
	rm -f ./funcs/*.html ./funcs/*.txt
