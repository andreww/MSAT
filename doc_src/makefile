default: html_funcs html_docs

AD=asciidoc 

asciidoc_funcs:
	(mkdir -p funcs; ./extract_matlab_help.py ../msat/*.m)

html_funcs: asciidoc_funcs
	(cd funcs ; for f in *.txt ; do $(AD) $$f ; done)

html_docs: MSATstart.txt MSATguide.txt
	$(AD) -a latexmath MSATstart.txt
	$(AD) MSATguide.txt

install: html_funcs html_docs
	mkdir -p ../msat/html/funcs
	cat helptoc.xml_header helptoc.xml_funcs helptoc.xml_footer > ../msat/html/helptoc.xml
	cp MSATstart.html MSATguide.html albite.png ../msat/html
	cp ./funcs/*.html ../msat/html/funcs
	$(MAKE) clean

external-install: html_funcs html_docs
	# NB - $(TARGET) must be set on the command line 
	# or this will eat your computer.
	mkdir -p $(TARGET)/funcs
	cp MSATstart.html MSATguide.html albite.png $(TARGET)
	cp ./funcs/*.html $(TARGET)/funcs
	$(MAKE) clean

clean:
	rm -f ./funcs/*.html ./funcs/*.txt helptoc.xml_funcs
	rm -f MSATstart.html MSATguide.html

clean_all: clean
	rm -rf ../msat/html/*	

