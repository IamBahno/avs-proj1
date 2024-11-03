set(CMAKE_C_COMPILER "/apps/all/intel-compilers/2024.2.0/compiler/2024.2/bin/icx")
set(CMAKE_C_COMPILER_ARG1 "")
set(CMAKE_C_COMPILER_ID "IntelLLVM")
set(CMAKE_C_COMPILER_VERSION "2024.2.0")
set(CMAKE_C_COMPILER_VERSION_INTERNAL "")
set(CMAKE_C_COMPILER_WRAPPER "")
set(CMAKE_C_STANDARD_COMPUTED_DEFAULT "17")
set(CMAKE_C_EXTENSIONS_COMPUTED_DEFAULT "ON")
set(CMAKE_C_COMPILE_FEATURES "c_std_90;c_function_prototypes;c_std_99;c_restrict;c_variadic_macros;c_std_11;c_static_assert;c_std_17")
set(CMAKE_C90_COMPILE_FEATURES "c_std_90;c_function_prototypes")
set(CMAKE_C99_COMPILE_FEATURES "c_std_99;c_restrict;c_variadic_macros")
set(CMAKE_C11_COMPILE_FEATURES "c_std_11;c_static_assert")
set(CMAKE_C17_COMPILE_FEATURES "c_std_17")
set(CMAKE_C23_COMPILE_FEATURES "")

set(CMAKE_C_PLATFORM_ID "Linux")
set(CMAKE_C_SIMULATE_ID "GNU")
set(CMAKE_C_COMPILER_FRONTEND_VARIANT "GNU")
set(CMAKE_C_SIMULATE_VERSION "4.2.1")




set(CMAKE_AR "/apps/all/binutils/2.42-GCCcore-13.3.0/bin/ar")
set(CMAKE_C_COMPILER_AR "")
set(CMAKE_RANLIB "/apps/all/binutils/2.42-GCCcore-13.3.0/bin/ranlib")
set(CMAKE_C_COMPILER_RANLIB "")
set(CMAKE_LINKER "/apps/all/binutils/2.42-GCCcore-13.3.0/bin/ld")
set(CMAKE_MT "")
set(CMAKE_COMPILER_IS_GNUCC )
set(CMAKE_C_COMPILER_LOADED 1)
set(CMAKE_C_COMPILER_WORKS TRUE)
set(CMAKE_C_ABI_COMPILED TRUE)

set(CMAKE_C_COMPILER_ENV_VAR "CC")

set(CMAKE_C_COMPILER_ID_RUN 1)
set(CMAKE_C_SOURCE_FILE_EXTENSIONS c;m)
set(CMAKE_C_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_C_LINKER_PREFERENCE 10)

# Save compiler ABI information.
set(CMAKE_C_SIZEOF_DATA_PTR "8")
set(CMAKE_C_COMPILER_ABI "ELF")
set(CMAKE_C_BYTE_ORDER "LITTLE_ENDIAN")
set(CMAKE_C_LIBRARY_ARCHITECTURE "x86_64-unknown-linux-gnu")

if(CMAKE_C_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_C_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_C_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_C_COMPILER_ABI}")
endif()

if(CMAKE_C_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "x86_64-unknown-linux-gnu")
endif()

set(CMAKE_C_CL_SHOWINCLUDES_PREFIX "")
if(CMAKE_C_CL_SHOWINCLUDES_PREFIX)
  set(CMAKE_CL_SHOWINCLUDES_PREFIX "${CMAKE_C_CL_SHOWINCLUDES_PREFIX}")
endif()





set(CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES "/apps/all/libarchive/3.6.1-GCCcore-11.3.0/include;/apps/all/XZ/5.2.5-GCCcore-11.3.0/include;/apps/all/cURL/7.83.0-GCCcore-11.3.0/include;/apps/all/OpenSSL/1.1/include;/apps/all/bzip2/1.0.8-GCCcore-11.3.0/include;/apps/all/zlib/1.2.12-GCCcore-11.3.0/include;/apps/all/ncurses/6.3-GCCcore-11.3.0/include;/apps/all/intel-compilers/2024.2.0/tbb/2021.13/include;/apps/all/binutils/2.42-GCCcore-13.3.0/include;/apps/all/Mesa/21.1.1-GCCcore-10.3.0/include;/apps/all/LLVM/11.1.0-GCCcore-10.3.0/include;/apps/all/libunwind/1.4.0-GCCcore-10.3.0/include;/apps/all/libglvnd/1.3.3-GCCcore-10.3.0/include;/apps/all/X11/20210518-GCCcore-10.3.0/include;/apps/all/libdrm/2.4.106-GCCcore-10.3.0/include;/apps/all/zstd/1.4.9-GCCcore-10.3.0/include;/apps/all/GTK2/2.24.33-GCCcore-10.3.0/include;/apps/all/Pango/1.48.5-GCCcore-10.3.0/include;/apps/all/FriBidi/1.0.10-GCCcore-10.3.0/include;/apps/all/HarfBuzz/2.8.1-GCCcore-10.3.0/include;/apps/all/cairo/1.16.0-GCCcore-10.3.0/include;/apps/all/freetype/2.10.4-GCCcore-10.3.0/include/freetype2;/apps/all/ICU/69.1-GCCcore-10.3.0/include;/apps/all/GLib/2.68.2-GCCcore-10.3.0/include;/apps/all/expat/2.2.9-GCCcore-10.3.0/include;/apps/all/pixman/0.40.0-GCCcore-10.3.0/include;/apps/all/libpng/1.6.37-GCCcore-10.3.0/include;/apps/all/Gdk-Pixbuf/2.42.6-GCCcore-10.3.0/include;/apps/all/libpciaccess/0.16-GCCcore-10.3.0/include;/apps/all/fontconfig/2.13.93-GCCcore-10.3.0/include;/apps/all/util-linux/2.36-GCCcore-10.3.0/include;/apps/all/Brotli/1.0.9-GCCcore-10.3.0/include;/apps/all/LibTIFF/4.2.0-GCCcore-10.3.0/include;/apps/all/lz4/1.9.3-GCCcore-10.3.0/include;/apps/all/jbigkit/2.1-GCCcore-10.3.0/include;/apps/all/libjpeg-turbo/2.0.6-GCCcore-10.3.0/include;/apps/all/ATK/2.36.0-GCCcore-10.3.0/include;/apps/all/PCRE2/10.36-GCCcore-10.3.0/include;/apps/all/gettext/0.21-GCCcore-10.3.0/include;/apps/all/libxml2/2.9.10-GCCcore-10.3.0/include/libxml2;/apps/all/libxml2/2.9.10-GCCcore-10.3.0/include;/apps/all/libffi/3.3-GCCcore-10.3.0/include;/apps/all/GTK3/3.24.31-GCCcore-11.2.0/include;/apps/all/libepoxy/1.5.8-GCCcore-11.2.0/include;/apps/all/at-spi2-atk/2.38.0-GCCcore-11.2.0/include;/apps/all/at-spi2-core/2.40.3-GCCcore-11.2.0/include;/apps/all/DBus/1.13.18-GCCcore-11.2.0/include;/apps/all/PCRE/8.45-GCCcore-11.2.0/include;/apps/all/imkl/2021.4.0/mkl/2021.4.0/include/fftw;/apps/all/imkl/2021.4.0/mkl/2021.4.0/include;/apps/all/impi/2021.4.0-intel-compilers-2021.4.0/mpi/2021.4.0/include;/apps/all/UCX/1.11.2-GCCcore-11.2.0/include;/apps/all/numactl/2.0.14-GCCcore-11.2.0/include;/apps/all/intel-compilers/2024.2.0/compiler/2024.2/opt/compiler/include;/apps/all/intel-compilers/2024.2.0/compiler/2024.2/lib/clang/19/include;/usr/local/include;/usr/include")
set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "svml;irng;imf;m;gcc;gcc_s;irc;dl;gcc;gcc_s;c;gcc;gcc_s;irc_s")
set(CMAKE_C_IMPLICIT_LINK_DIRECTORIES "/apps/all/intel-compilers/2024.2.0/compiler/2024.2/lib;/apps/all/intel-compilers/2024.2.0/compiler/2024.2/lib/clang/19/lib/x86_64-unknown-linux-gnu;/apps/all/GCCcore/11.3.0/lib/gcc/x86_64-pc-linux-gnu/11.3.0;/apps/all/GCCcore/11.3.0/lib64;/lib64;/usr/lib64;/apps/all/GCCcore/11.3.0/lib;/apps/all/intel-compilers/2024.2.0/compiler/2024.2/opt/compiler/lib;/lib;/usr/lib;/apps/all/libarchive/3.6.1-GCCcore-11.3.0/lib;/apps/all/XZ/5.2.5-GCCcore-11.3.0/lib;/apps/all/cURL/7.83.0-GCCcore-11.3.0/lib;/apps/all/OpenSSL/1.1/lib;/apps/all/bzip2/1.0.8-GCCcore-11.3.0/lib;/apps/all/zlib/1.2.12-GCCcore-11.3.0/lib;/apps/all/ncurses/6.3-GCCcore-11.3.0/lib;/apps/all/intel-compilers/2024.2.0/tbb/2021.13/lib/intel64/gcc4.8;/apps/all/binutils/2.42-GCCcore-13.3.0/lib;/apps/all/Mesa/21.1.1-GCCcore-10.3.0/lib;/apps/all/LLVM/11.1.0-GCCcore-10.3.0/lib;/apps/all/libunwind/1.4.0-GCCcore-10.3.0/lib;/apps/all/libglvnd/1.3.3-GCCcore-10.3.0/lib;/apps/all/X11/20210518-GCCcore-10.3.0/lib;/apps/all/libdrm/2.4.106-GCCcore-10.3.0/lib;/apps/all/zstd/1.4.9-GCCcore-10.3.0/lib;/apps/all/GTK2/2.24.33-GCCcore-10.3.0/lib;/apps/all/Pango/1.48.5-GCCcore-10.3.0/lib;/apps/all/FriBidi/1.0.10-GCCcore-10.3.0/lib;/apps/all/HarfBuzz/2.8.1-GCCcore-10.3.0/lib;/apps/all/cairo/1.16.0-GCCcore-10.3.0/lib;/apps/all/freetype/2.10.4-GCCcore-10.3.0/lib;/apps/all/ICU/69.1-GCCcore-10.3.0/lib;/apps/all/GLib/2.68.2-GCCcore-10.3.0/lib;/apps/all/expat/2.2.9-GCCcore-10.3.0/lib;/apps/all/pixman/0.40.0-GCCcore-10.3.0/lib;/apps/all/libpng/1.6.37-GCCcore-10.3.0/lib;/apps/all/Gdk-Pixbuf/2.42.6-GCCcore-10.3.0/lib;/apps/all/libpciaccess/0.16-GCCcore-10.3.0/lib;/apps/all/fontconfig/2.13.93-GCCcore-10.3.0/lib;/apps/all/util-linux/2.36-GCCcore-10.3.0/lib;/apps/all/Brotli/1.0.9-GCCcore-10.3.0/lib;/apps/all/LibTIFF/4.2.0-GCCcore-10.3.0/lib;/apps/all/lz4/1.9.3-GCCcore-10.3.0/lib;/apps/all/jbigkit/2.1-GCCcore-10.3.0/lib;/apps/all/libjpeg-turbo/2.0.6-GCCcore-10.3.0/lib;/apps/all/ATK/2.36.0-GCCcore-10.3.0/lib;/apps/all/PCRE2/10.36-GCCcore-10.3.0/lib;/apps/all/gettext/0.21-GCCcore-10.3.0/lib;/apps/all/libxml2/2.9.10-GCCcore-10.3.0/lib;/apps/all/libffi/3.3-GCCcore-10.3.0/lib64;/apps/all/libffi/3.3-GCCcore-10.3.0/lib;/apps/all/GTK3/3.24.31-GCCcore-11.2.0/lib;/apps/all/libepoxy/1.5.8-GCCcore-11.2.0/lib;/apps/all/at-spi2-atk/2.38.0-GCCcore-11.2.0/lib;/apps/all/at-spi2-core/2.40.3-GCCcore-11.2.0/lib;/apps/all/DBus/1.13.18-GCCcore-11.2.0/lib;/apps/all/PCRE/8.45-GCCcore-11.2.0/lib;/apps/all/imkl-FFTW/2021.4.0-iimpi-2021b/lib;/apps/all/imkl/2021.4.0/mkl/2021.4.0/lib/intel64;/apps/all/imkl/2021.4.0/compiler/2021.4.0/linux/compiler/lib/intel64_lin;/apps/all/impi/2021.4.0-intel-compilers-2021.4.0/mpi/2021.4.0/libfabric/lib;/apps/all/impi/2021.4.0-intel-compilers-2021.4.0/mpi/2021.4.0/lib/release;/apps/all/impi/2021.4.0-intel-compilers-2021.4.0/mpi/2021.4.0/lib;/apps/all/UCX/1.11.2-GCCcore-11.2.0/lib;/apps/all/numactl/2.0.14-GCCcore-11.2.0/lib")
set(CMAKE_C_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
