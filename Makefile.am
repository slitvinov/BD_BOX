SUBDIRS = src
dist_doc_DATA = README
EXTRA_DIST=examples src/win \
           src/bd_box.mak src/bd_box.sln \
           src/bd_box_gpu.mak src/bd_box_gpu.sln \
           src/erfc.c
dist-hook:
	rm -rf `find $(distdir)/examples -name .svn`
	rm -rf `find $(distdir)/src/win -name .svn`
