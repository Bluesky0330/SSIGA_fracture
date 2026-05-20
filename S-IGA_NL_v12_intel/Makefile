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

# 四点曲げ試験/溶接協会, small/finite, iga/siga, coarse/fine ep
4tenmage_fine_iga: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_IGA_ep_katayama.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/pipe_input.txt ./option_IGA_ep_katayama.ini $(COMMANDS2)

4tenmage_fine_siga: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_SIGA_ep_katayama.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/pipe_input_glo.txt ../input/pipe_input_loc.txt ./option_SIGA_ep_katayama.ini $(COMMANDS2)

# omar
omar_1: compile rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_crack_J_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src/$(PROGRAM_NAME).x ../input/omar_input_Embedded_crack_a=2_c=2_singular_crack.txt ./option_crack_J_e.ini $(COMMANDS2)
# 	cd ./analysis/$(OUTPUT_DIR) && valgrind --tool=memcheck --leak-check=full ../../src/$(PROGRAM_NAME).x ../input/omar_input_Embedded_crack_a=2_c=2_singular_crack.txt ./option_crack_J_e.ini

omar_1_optimized: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_crack_J_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/omar_input_Embedded_crack_a=2_c=2_singular_crack.txt ./option_crack_J_e.ini $(COMMANDS2)
# 	cd ./analysis/$(OUTPUT_DIR) && valgrind --tool=memcheck --leak-check=full ../../src/$(PROGRAM_NAME).x ../input/omar_input_Embedded_crack_a=2_c=2_singular_crack.txt ./option_crack_J_e.ini

omar_2: compile rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_crack_J_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src/$(PROGRAM_NAME).x ../input/omar_input_Embedded_crack_a=1.2_c=2_singular_crack.txt ./option_crack_J_e.ini $(COMMANDS2)

omar_3: compile rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_crack_J_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src/$(PROGRAM_NAME).x ../input/omar_input_Embedded_crack_a=0.4_c=2_singular_crack_Ver5.txt ./option_crack_J_e.ini $(COMMANDS2)

omar_ep_1: compile rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_crack_J_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src/$(PROGRAM_NAME).x ../input/omar_input_Embedded_crack_a=2_c=2_singular_crack.txt ./option_crack_J_ep.ini $(COMMANDS2)

omar_ep_2: compile rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_crack_J_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src/$(PROGRAM_NAME).x ../input/omar_input_Embedded_crack_a=1.2_c=2_singular_crack.txt ./option_crack_J_ep.ini $(COMMANDS2)

omar_ep_3: compile rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_crack_J_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src/$(PROGRAM_NAME).x ../input/omar_input_Embedded_crack_a=0.4_c=2_singular_crack_Ver5.txt ./option_crack_J_ep.ini $(COMMANDS2)

# sample
IGA_through_crack_e: compile rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_IGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src/$(PROGRAM_NAME).x ../input/IGA_through_crack.txt ./option_IGA_e.ini $(COMMANDS2)

IGA_through_crack_e_optimized: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_IGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/IGA_through_crack.txt ./option_IGA_e.ini $(COMMANDS2)

IGA_through_crack_ep: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_IGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/IGA_through_crack.txt ./option_IGA_ep.ini $(COMMANDS2)

IGA_through_crack_4th_order_ep: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_IGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/IGA_through_crack_4th_order.txt ./option_IGA_ep.ini $(COMMANDS2)

IGA_through_crack_5th_order_ep: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_IGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/IGA_through_crack_5th_order.txt ./option_IGA_ep.ini $(COMMANDS2)

IGA_SCCHP: compile rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_IGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src/$(PROGRAM_NAME).x ../input/SCCHP.txt ./option_IGA_e.ini $(COMMANDS2)

# APCOM
APCOM_IGA_coarse1_1.0_e: compile rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_IGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src/$(PROGRAM_NAME).x ../input/APCOM_iga_coarse1_1.0.txt ./option_APCOM_IGA_e.ini $(COMMANDS2)

APCOM_SIGA_coarse1_1.0_e: compile rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src/$(PROGRAM_NAME).x ../input/APCOM_glo_coarse1_1.0.txt ../input/APCOM_loc_coarse1_1.0.txt ./option_APCOM_SIGA_e.ini $(COMMANDS2)

APCOM_SIGA_coarse1_1.0_e_op: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_coarse1_1.0.txt ../input/APCOM_loc_coarse1_1.0.txt ./option_APCOM_SIGA_e.ini $(COMMANDS2)

APCOM_SIGA_coarse1_1.0_ep: compile rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src/$(PROGRAM_NAME).x ../input/APCOM_glo_coarse1_1.0.txt ../input/APCOM_loc_coarse1_1.0.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_IGA_fine1_1.0_e: compile rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_IGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src/$(PROGRAM_NAME).x ../input/APCOM_iga_fine1_1.0.txt ./option_APCOM_IGA_e.ini $(COMMANDS2)

APCOM_IGA_fine1_1.0_ep: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_IGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_iga_fine1_1.0.txt ./option_APCOM_IGA_ep.ini $(COMMANDS2)

APCOM_IGA_fine1_0.8_ep: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_IGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_iga_fine1_0.8.txt ./option_APCOM_IGA_ep.ini $(COMMANDS2)

APCOM_IGA_fine1_0.6_ep: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_IGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_iga_fine1_0.6.txt ./option_APCOM_IGA_ep.ini $(COMMANDS2)

APCOM_IGA_fine1_0.4_ep: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_IGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_iga_fine1_0.4.txt ./option_APCOM_IGA_ep.ini $(COMMANDS2)

APCOM_IGA_fine1_0.2_ep: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_IGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_iga_fine1_0.2.txt ./option_APCOM_IGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_1.0_ep: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_1.0.txt ../input/APCOM_loc_fine1_1.0.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.8_ep: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.8.txt ../input/APCOM_loc_fine1_0.8.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.6_ep: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.6.txt ../input/APCOM_loc_fine1_0.6.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.4_ep: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.4.txt ../input/APCOM_loc_fine1_0.4.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.2_ep: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.2.txt ../input/APCOM_loc_fine1_0.2.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_IGA_large_deformation: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_IGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_iga_large_deformation.txt ./option_APCOM_IGA_ep.ini $(COMMANDS2)

APCOM_SIGA_large_deformation: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_large_deformation.txt ../input/APCOM_loc_large_deformation.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_IGA_large_deformation_through_crack: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_TC_IGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_iga_large_deformation_through_crack.txt ./option_APCOM_TC_IGA_ep.ini $(COMMANDS2)

APCOM_SIGA_large_deformation_through_crack: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_TC_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_large_deformation_through_crack.txt ../input/APCOM_loc_large_deformation_through_crack.txt ./option_APCOM_TC_SIGA_ep.ini $(COMMANDS2)

hamano_SIGA_e: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_SIGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/hamano_glo.txt ../input/hamano_loc.txt ./option_SIGA_e.ini $(COMMANDS2)

hamano_coarse_SIGA_e: compile rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_SIGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src/$(PROGRAM_NAME).x ../input/hamano_coarse_glo.txt ../input/hamano_coarse_loc.txt ./option_SIGA_e.ini $(COMMANDS2)

simple_crack: compile rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_SIGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src/$(PROGRAM_NAME).x ../input/simple_crack_glo.txt ../input/simple_crack_loc.txt ./option_SIGA_e.ini $(COMMANDS2)

singular_patch_test: compile rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_SIGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src/$(PROGRAM_NAME).x ../input/singular_patch_test_glo.txt ../input/singular_patch_test_loc.txt ./option_SIGA_e.ini $(COMMANDS2)

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