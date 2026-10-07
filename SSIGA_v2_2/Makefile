PROGRAM_NAME = SSIGA

# INPUT_NAME := pwh_IGA
# INPUT_NAME := pwh_loc_8x8
# INPUT_NAME := pwh_loc_16x16
# INPUT_NAME := pwh_loc_16x16_OG
# INPUT_NAME := pwh_loc_32x32
# INPUT_NAME := pwh_loc_32x32_OG
# INPUT_NAME := SCCHP_v3_coarse
# INPUT_NAME := SCCHP_v3_fine_OG
# INPUT_NAME := SCCHP_v3_fine
# INPUT_NAME := SCCHP_v6_coarse
# INPUT_NAME := SCCHP_v6_fine
# INPUT_NAME := SCCHP_v10_coarse
# INPUT_NAME := SCCHP_v10_fine
# INPUT_NAME := SCCHP_v11_coarse
# INPUT_NAME := SCCHP_v11_fine
# INPUT_NAME := SCCHP_v12_coarse
# INPUT_NAME := SCCHP_v12_fine
# INPUT_NAME := SCCHP_v13_coarse
# INPUT_NAME := SCCHP_v13_fine
# INPUT_NAME := SCCHP_v14_coarse
# INPUT_NAME := SCCHP_v14_fine
# INPUT_NAME := SCCHP_v15_coarse
# INPUT_NAME := SCCHP_v15_fine
# INPUT_NAME := SCCHP_v16_coarse
# INPUT_NAME := SCCHP_v16_fine
# INPUT_NAME := SCCHP_v17_coarse
# INPUT_NAME := SCCHP_v17_fine
# INPUT_NAME := SCCHP_v18_coarse
# INPUT_NAME := SCCHP_v18_fine
# INPUT_NAME := SCCHP_v19_coarse
# INPUT_NAME := SCCHP_v19_fine
# INPUT_NAME := PBTC_IGA
# INPUT_NAME := PBTC_test
# INPUT_NAME := PBTC_v0_coarse
# INPUT_NAME := PBTC_v0_fine
# INPUT_NAME := PBTC_v1_coarse
# INPUT_NAME := PBTC_v1_fine
# INPUT_NAME := PBTC_v1_sf
# INPUT_NAME := PBTC_v1_sc
# INPUT_NAME := PBTC_v1_fine_OG
# INPUT_NAME := PBTC_v1_order_diff
# INPUT_NAME := EPBTC_IGA
# INPUT_NAME := beam

# 出力ディレクトリをインプット名から自動生成
OUTPUT_DIR := output_$(INPUT_NAME)

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

SSIGA: compile rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_SSIGA.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src/$(PROGRAM_NAME).x ../input/$(INPUT_NAME)/input_glo.txt ../input/$(INPUT_NAME)/input_loc.txt ./option_SSIGA.ini $(COMMANDS2)

SSIGA_optimized: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_SSIGA.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/$(INPUT_NAME)/input_glo.txt ../input/$(INPUT_NAME)/input_loc.txt ./option_SSIGA.ini $(COMMANDS2)

IGA: compile rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_IGA.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src/$(PROGRAM_NAME).x ../input/$(INPUT_NAME)/input.txt ./option_IGA.ini $(COMMANDS2)

debug: compile rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_SSIGA.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src/$(PROGRAM_NAME).x ../input/$(INPUT_NAME)/input_glo.txt ../input/$(INPUT_NAME)/input_loc.txt ./option_SSIGA.ini $(COMMANDS2)

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
	@echo "Cleaning $(OUTPUT_DIR) while preserving __Displacement.dat"
	@cd ./analysis && \
	TMP_FILE=".keep___Displacement.dat" && \
	if [ -f "$(OUTPUT_DIR)/__Displacement.dat" ]; then cp -f "$(OUTPUT_DIR)/__Displacement.dat" "$$TMP_FILE"; echo "Preserved __Displacement.dat"; fi && \
	rm -rf "$(OUTPUT_DIR)" && \
	mkdir -p "$(OUTPUT_DIR)" && \
	if [ -f "$$TMP_FILE" ]; then mv -f "$$TMP_FILE" "$(OUTPUT_DIR)/__Displacement.dat"; fi

# clean
clean: rm
	rm -f ./src/*.x
	rm -f ./src/*.o
	rm -f ./src/*.a
	rm -f ./src/*.d
	rm -f ./src/*.tmp