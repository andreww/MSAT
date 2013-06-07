default: html_funcs html_docs

AD=asciidoc 

asciidoc_funcs:
	(mkdir -p funcs; ./extract_matlab_help.py ../msat/*.m)

html_funcs: asciidoc_funcs
	(mkdir -p funcs; cd funcs ; for f in *.txt ; do $(AD) $$f ; done)

html_docs: MSATstart.txt MSATguide.txt MSATfunctions.txt ../RELEASE_NOTES example-docs
	$(AD) -a latexmath MSATstart.txt
	$(AD) MSATguide.txt
	$(AD) MSATdevelop.txt
	$(AD) -a latexmath MSATfunctions.txt
	$(AD) -o MSATrelease_notes.html ../RELEASE_NOTES 
	$(AD) -o MSATreadme.html ../README

install: html_funcs html_docs example-docs
	mkdir -p ../msat/html/funcs
	mkdir -p ../msat/html/examples
	cat helptoc.xml_header helptoc.xml_funcs helptoc.xml_footer > ../msat/html/helptoc.xml
	cp MSATdevelop.html MSATrelease_notes.html MSATreadme.html MSATfunctions.html  MSATstart.html MSATguide.html albite.png MSATfunctions.html ../msat/html
	cp ./funcs/*.html ../msat/html/funcs
	cp ./examples/*.html ../msat/html/examples
	$(MAKE) clean

example-docs: examples/diopside_example.txt examples/index.txt examples/lmm_anisotropy_example.txt examples/loading_examples.txt examples/splitting_model_example.txt examples/texture_example.txt
	(cd examples; $(AD) diopside_example.txt; $(AD) index.txt ; $(AD) lmm_anisotropy_example.txt ; $(AD) loading_examples.txt ; $(AD) splitting_model_example.txt ; $(AD) texture_example.txt)

external-install: html_funcs html_docs
	# NB - $(TARGET) must be set on the command line 
	# or this will eat your computer.
	mkdir -p $(TARGET)/funcs
	mkdir -p $(TARGET)/examples
	cp MSATstart.html MSATdevelop.html MSATguide.html MSATreadme.html MSATrelease_notes.html MSATfunctions.html albite.png $(TARGET)
	cp ./funcs/*.html $(TARGET)/funcs
	cp ./examples/*.html $(TARGET)/examples
	$(MAKE) clean

clean:
	rm -f ./funcs/*.html ./funcs/*.txt helptoc.xml_funcs
	rm -f ./examples/*.html
	rm -f MSATstart.html MSATguide.html MSATfunctions.html
	rm -f MSATrelease_notes.html MSATreadme.html MSATdevelop.html
	rm -rf funcs

clean_all: clean
	rm -rf ../msat/html/*	

