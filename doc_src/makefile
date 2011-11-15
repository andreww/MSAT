default: html_funcs html_docs

asciidoc_funcs:
	./extract_matlab_help.py ../msat/*.m

html_funcs: asciidoc_funcs
	(cd funcs ; for f in *.txt ; do asciidoc $$f ; done)

html_docs: MSATstart.txt MSATguide.txt
	asciidoc MSATstart.txt
	asciidoc MSATguide.txt

install: html_funcs html_docs
	mkdir -p ../msat/html/funcs
	cat helptoc.xml_header helptoc.xml_funcs helptoc.xml_footer > ../msat/html/helptoc.xml
	cp MSATstart.html MSATguide.html ../msat/html
	cp ./funcs/*.html ../msat/html/funcs

clean:
	rm -f ./funcs/*.html ./funcs/*.txt helptoc.xml_funcs
	rm -f MSATstart.html MSATguide.html

clean_all: clean
	rm -rf ../msat/html/*	

