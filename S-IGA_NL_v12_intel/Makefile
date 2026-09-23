PROGRAM_NAME	= SIGA

# 出力ディレクトリをターゲット名から自動生成
OUTPUT_DIR := output_$(MAKECMDGOALS)

NUMACTL := $(shell numactl --hardware > /dev/null 2>&1 && echo numactl --localalloc || echo "")
PERF_STAT := $(shell perf --version > /dev/null 2>&1 && echo perf stat -dd || echo "")
PERF_RECORD := $(shell perf --version > /dev/null 2>&1 && echo perf record || echo "")
STDBUF := $(shell which stdbuf > /dev/null 2>&1 && echo stdbuf -oL || echo "")
TEE := $(shell which tee > /dev/null 2>&1 && echo tee || echo "")

COMMANDS1 :=
COMMANDS2 :=

override PERF_MODE = 0 # 0: stat, 1: record

PERF_MODE_STRIPPED := $(strip $(PERF_MODE))
ifeq ($(PERF_MODE_STRIPPED),0)
    $(info Perf mode: stat)
    COMMANDS1 := $(PERF_STAT) $(NUMACTL) $(STDBUF)
    COMMANDS2 := | $(TEE) log.txt
else
    $(info Perf mode: record)
    COMMANDS1 := $(PERF_RECORD) $(NUMACTL)
    COMMANDS2 := && perf report --stdio -i ./analysis/$(OUTPUT_DIR)/perf.data > ./analysis/$(OUTPUT_DIR)/perf_report.txt && head -25 ./analysis/$(OUTPUT_DIR)/perf_report.txt
endif

# compile
compile:
	cd ./src && make -f Make_src

# compile
compile_optimized:
	cd ./src_fast && ./build.sh

# plate with hole
pwh: compile rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_SIGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src/$(PROGRAM_NAME).x ../input/pwh_glo_8x8.txt ../input/pwh_loc_16x16.txt ./option_SIGA_e.ini $(COMMANDS2)

# for paraview
view_0:
	cd ./analysis/$(OUTPUT_DIR) && paraview global_patch.xmf2 global_patch_boundary_line.xmf2 global_patch_control_point.xmf2

view_1:
	cd ./analysis/$(OUTPUT_DIR) && paraview global_patch.xmf2 global_patch_boundary_line.xmf2 global_patch_control_point.xmf2 \
	local_patch.xmf2 local_patch_boundary_line.xmf2 local_patch_control_point.xmf2

view_2:
	cd ./analysis/$(OUTPUT_DIR) && paraview global_patch.xmf2 global_patch_boundary_line.xmf2 global_patch_control_point.xmf2 \
	deformed_global_patch.xmf2 deformed_global_patch_boundary_line.xmf2

view_3:
	cd ./analysis/$(OUTPUT_DIR) && paraview global_patch.xmf2 global_patch_boundary_line.xmf2 global_patch_control_point.xmf2 \
	local_patch.xmf2 local_patch_boundary_line.xmf2 local_patch_control_point.xmf2 \
	deformed_global_patch.xmf2 deformed_global_patch_boundary_line.xmf2 \
	deformed_local_patch.xmf2 deformed_local_patch_boundary_line.xmf2

# mkdir
mkdir_output:
	cd ./analysis && mkdir -p $(OUTPUT_DIR)
	cd ./analysis/$(OUTPUT_DIR) && mkdir -p bin

# rm
rm:
	cd ./analysis && rm -rf $(OUTPUT_DIR) && mkdir -p $(OUTPUT_DIR)/bin

# clean
clean: rm
	rm -f ./src/*.x
	rm -f ./src/*.o
	rm -f ./src/*.a
	rm -f ./src/*.d
	rm -f ./src/*.tmp
	rm -f ./src_fast/src/*.x
	rm -f ./src_fast/src/*.o
	rm -f ./src_fast/src/*.a
	rm -f ./src_fast/src/*.d
	rm -f ./src_fast/src/*.tmp
	rm -f ./src_fast/*.x