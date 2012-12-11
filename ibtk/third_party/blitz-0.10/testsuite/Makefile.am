#
# Written by Patrick Guio <patrick.guio@fys.uio.no>
#

# Blitz++ test suite makefile

EXTRA_DIST = testsuite.h

AM_CPPFLAGS = -I$(srcdir) -I$(top_srcdir) -I$(top_builddir) $(BOOST_CPPFLAGS)
AM_CXXFLAGS = @CXX_DEBUG_FLAGS@ -DBZ_DEBUG
AM_LDFLAGS = $(BOOST_LDFLAGS)
LDADD = -L$(top_builddir)/lib -lblitz $(BOOST_SERIALIZATION_LIB)

EXTRA_PROGRAMS = 64bit Adnene-Ben-Abdallah-1 Adnene-Ben-Abdallah-2     \
alignment arrayresize arrayinitialize bitwise chris-jeffery-1 chris-jeffery-2  \
chris-jeffery-3 complex-test constarray contiguous copy ctors	       \
derrick-bass-1 derrick-bass-3 exprctor expression-slicing extract free \
gary-huber-1 indexexpr-base initialize int-math-func interlace iter    \
Josef-Wagenhuber levicivita loop1 matthias-troyer-1 matthias-troyer-2  \
mattias-lindstroem-1 member_function minmax minsumpow module	       \
multicomponent multicomponent-2 newet Olaf-Ronneberger-1	       \
patrik-jonsson-1 peter-bienstman-1 peter-bienstman-2 peter-bienstman-3 \
peter-bienstman-4 peter-bienstman-5 peter-nordlund-1 peter-nordlund-2  \
peter-nordlund-3 preexisting promote pthread qcd reduce reindex	       \
reverse safeToReturn shapecheck shape slice-iterators stencil-et       \
stencil-extent storage stub theodore-papadopoulo-1 tinymat tinyvec     \
transpose troyer-genilloud tvinitialize Ulisses-Mello-1 weakref	       \
wei-ku-1 where zeek-1 cast serialize update

# 64bit test is special because it takes forever to run without optimization
64bit_SOURCES = 64bit.cpp
64bit_CXXFLAGS = @CXX_OPTIMIZE_FLAGS@ -DBZ_DEBUG

Adnene_Ben_Abdallah_1_SOURCES = Adnene-Ben-Abdallah-1.cpp
Adnene_Ben_Abdallah_2_SOURCES = Adnene-Ben-Abdallah-2.cpp
alignment_SOURCES = alignment.cpp
arrayinitialize_SOURCES = arrayinitialize.cpp
arrayresize_SOURCES = arrayresize.cpp
bitwise_SOURCES = bitwise.cpp
cast_SOURCES = cast.cpp
chris_jeffery_1_SOURCES = chris-jeffery-1.cpp
chris_jeffery_2_SOURCES = chris-jeffery-2.cpp
chris_jeffery_3_SOURCES = chris-jeffery-3.cpp
complex_test_SOURCES = complex-test.cpp
constarray_SOURCES = constarray.cpp
contiguous_SOURCES = contiguous.cpp
copy_SOURCES = copy.cpp
ctors_SOURCES = ctors.cpp
derrick_bass_1_SOURCES = derrick-bass-1.cpp
derrick_bass_3_SOURCES = derrick-bass-3.cpp
exprctor_SOURCES = exprctor.cpp
expression_slicing_SOURCES = expression-slicing.cpp
extract_SOURCES = extract.cpp
free_SOURCES = free.cpp
gary_huber_1_SOURCES = gary-huber-1.cpp
indexexpr_base_SOURCES = indexexpr-base.cpp
initialize_SOURCES = initialize.cpp
int_math_func_SOURCES = int-math-func.cpp
interlace_SOURCES = interlace.cpp
iter_SOURCES = iter.cpp
Josef_Wagenhuber_SOURCES = Josef-Wagenhuber.cpp
levicivita_SOURCES = levicivita.cpp
loop1_SOURCES = loop1.cpp
matthias_troyer_1_SOURCES = matthias-troyer-1.cpp
matthias_troyer_2_SOURCES = matthias-troyer-2.cpp
mattias_lindstroem_1_SOURCES = mattias-lindstroem-1.cpp
member_function_SOURCES = member_function.cpp
minmax_SOURCES = minmax.cpp
minsumpow_SOURCES = minsumpow.cpp
module_SOURCES = module1.cpp module2.cpp
multicomponent_SOURCES = multicomponent.cpp
multicomponent_2_SOURCES = multicomponent-2.cpp
newet_SOURCES = newet.cpp
Olaf_Ronneberger_1_SOURCES = Olaf-Ronneberger-1.cpp
patrik_jonsson_1_SOURCES = patrik-jonsson-1.cpp
peter_bienstman_1_SOURCES = peter-bienstman-1.cpp
peter_bienstman_2_SOURCES = peter-bienstman-2.cpp
peter_bienstman_3_SOURCES = peter-bienstman-3.cpp
peter_bienstman_4_SOURCES = peter-bienstman-4.cpp
peter_bienstman_5_SOURCES = peter-bienstman-5.cpp
peter_nordlund_1_SOURCES = peter-nordlund-1.cpp
peter_nordlund_2_SOURCES = peter-nordlund-2.cpp
peter_nordlund_3_SOURCES = peter-nordlund-3.cpp
preexisting_SOURCES = preexisting.cpp
promote_SOURCES = promote.cpp
pthread_SOURCES = pthread.cpp
qcd_SOURCES = qcd.cpp
reduce_SOURCES = reduce.cpp
reindex_SOURCES = reindex.cpp
reverse_SOURCES = reverse.cpp
safeToReturn_SOURCES = safeToReturn.cpp
serialize_SOURCES = serialize.cpp
shapecheck_SOURCES = shapecheck.cpp
shape_SOURCES = shape.cpp
slice_iterators_SOURCES = slice-iterators.cpp
stencil_et_SOURCES = stencil-et.cpp
stencil_extent_SOURCES = stencil-extent.cpp
storage_SOURCES = storage.cpp
stub_SOURCES = stub.cpp
theodore_papadopoulo_1_SOURCES = theodore-papadopoulo-1.cpp
tinymat_SOURCES = tinymat.cpp
tinyvec_SOURCES = tinyvec.cpp
transpose_SOURCES = transpose.cpp
troyer_genilloud_SOURCES = troyer-genilloud.cpp
tvinitialize_SOURCES = tvinitialize.cpp
Ulisses_Mello_1_SOURCES = Ulisses-Mello-1.cpp
update_SOURCES = update.cpp
weakref_SOURCES = weakref.cpp
wei_ku_1_SOURCES = wei-ku-1.cpp
where_SOURCES = where.cpp
zeek_1_SOURCES = zeek-1.cpp

check-testsuite:  $(EXTRA_PROGRAMS)
	@echo Running test suite...
	@fail=0; for prog in $(EXTRA_PROGRAMS) ; do \
		if (./$$prog > /dev/null 2>&1); then \
			echo $$prog... passed; \
		else \
			echo TEST FAILED: $$prog; fail=1; fi; \
	done; \
	if [[ $$fail -eq 1 ]]; then echo There were failing tests; else echo All tests passed; fi

if HAVE_VALGRIND
valgrind-testsuite:  $(EXTRA_PROGRAMS)
	@echo Running test suite through valgrind...
	@for prog in $(EXTRA_PROGRAMS) ; do \
		echo -en "$$prog\t" && (valgrind --log-file-exactly=$$prog.val ./$$prog > /dev/null 2>&1 ;\
                grep "ERROR SUMMARY" $$prog.val) || exit 1; \
	done
	@echo All tests run.
endif

clean-local:
	-rm -rf *.ii *.ti cxx_repository Template.dir ii_files ti_files $(EXTRA_PROGRAMS)

