# ----------------------------------------------------
# CPUフラグの検出
# ----------------------------------------------------
CPU_FLAGS := $(shell grep -m1 -o -E 'avx512f|avx2|avx' /proc/cpuinfo | sort -u)
ARCH_OPT  :=

# ----------------------------------------------------
# オプションの分岐設定
# ----------------------------------------------------
ifeq ($(findstring avx512f,$(CPU_FLAGS)),avx512f)
    ARCH_OPT := -qopt-zmm-usage=high
else ifeq ($(findstring avx2,$(CPU_FLAGS)),avx2)
    ARCH_OPT := -qopt-zmm-usage=low
else
    ARCH_OPT :=
endif
# ----------------------------------------------------

PROGRAM_NAME = SIGA
COMPILER     = icpx

# ----------------------------------------------------
# 最適化オプション設定
# ----------------------------------------------------
COMMON_OPT := -O3 -march=native -mtune=native -fp-model precise \
              -qopenmp -qopenmp-simd \
              -qopt-prefetch \
              -qopt-mem-layout-trans=3 \
              -fno-alias \
              -falign-functions=64 -std=c++20 \
              -fdiagnostics-color=always \
              $(ARCH_OPT)

# MKLの設定（並列版）
MKL_OPT    := -qmkl=parallel

# LTO（Link Time Optimization）の設定
LTO_OPT    := -flto=auto -flto-jobs=$(shell nproc)

merged_output.cpp: constant.ini main.cpp class.cpp sub_MI3D.cpp sub_post.cpp sub_pre.cpp sub_search.cpp sub_solver.cpp sub_tool.cpp
	../fast.x ./constant.ini 3 ./main.cpp ./class.cpp ./sub_MI3D.cpp ./sub_post.cpp ./sub_pre.cpp ./sub_search.cpp ./sub_solver.cpp ./sub_tool.cpp

# ----------------------------------------------------
# コンパイル・リンク
# ----------------------------------------------------
$(PROGRAM_NAME).x: merged_output.cpp
	$(COMPILER) $(COMMON_OPT) $(MKL_OPT) $(LTO_OPT) -o $(PROGRAM_NAME).x merged_output.cpp \
        -ipo \
        -Wl,--start-group \
        -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core \
        -Wl,--end-group \
        -liomp5 -lpthread -lm -ldl \
        -lhdf5_cpp -lhdf5 -L/usr/lib/x86_64-linux-gnu/hdf5/serial \
        -Wall -Wextra -I/usr/include/hdf5/serial