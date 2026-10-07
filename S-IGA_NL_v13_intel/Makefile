PROGRAM_NAME	= SIGA

# 出力ディレクトリをターゲット名から自動生成
OUTPUT_DIR := output_$(MAKECMDGOALS)

NUMACTL := $(shell numactl --hardware > /dev/null 2>&1 && echo numactl --localalloc || echo "")
PERF_STAT := $(shell perf --version > /dev/null 2>&1 && echo perf stat -dd || echo "")
PERF_RECORD := $(shell perf --version > /dev/null 2>&1 && echo perf record || echo "")
STDBUF := $(shell which stdbuf > /dev/null 2>&1 && echo stdbuf -oL || echo "")
TEE := $(shell which tee > /dev/null 2>&1 && echo tee || echo "")
TIME_CMD := $(shell which time > /dev/null 2>&1 && echo /usr/bin/time -v || echo "")

COMMANDS1 :=
COMMANDS2 :=

override PERF_MODE = 0 # 0: stat, 1: record

PERF_MODE_STRIPPED := $(strip $(PERF_MODE))
ifeq ($(PERF_MODE_STRIPPED),0)
    $(info Perf mode: stat)
    COMMANDS1 := $(PERF_STAT) $(TIME_CMD) $(NUMACTL) $(STDBUF)
    COMMANDS2 := | $(TEE) log.txt
else
    $(info Perf mode: record)
    COMMANDS1 := $(PERF_RECORD) $(TIME_CMD) $(NUMACTL)
    COMMANDS2 := && perf report --stdio -i ./analysis/$(OUTPUT_DIR)/perf.data > ./analysis/$(OUTPUT_DIR)/perf_report.txt && head -25 ./analysis/$(OUTPUT_DIR)/perf_report.txt
endif

# compile
compile:
	cd ./src && make -j$(nproc) -f Make_src

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

APCOM_IGA_fine1_1.0_e: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_IGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_iga_fine1_1.0.txt ./option_APCOM_IGA_e.ini $(COMMANDS2)

APCOM_IGA_fine1_0.8_e: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_IGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_iga_fine1_0.8.txt ./option_APCOM_IGA_e.ini $(COMMANDS2)

APCOM_IGA_fine1_0.6_e: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_IGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_iga_fine1_0.6.txt ./option_APCOM_IGA_e.ini $(COMMANDS2)

APCOM_IGA_fine1_0.4_e: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_IGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_iga_fine1_0.4.txt ./option_APCOM_IGA_e.ini $(COMMANDS2)

APCOM_IGA_fine1_0.2_e: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_IGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_iga_fine1_0.2.txt ./option_APCOM_IGA_e.ini $(COMMANDS2)

APCOM_SIGA_fine1_1.0_e: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_1.0.txt ../input/APCOM_loc_fine1_1.0.txt ./option_APCOM_SIGA_e.ini $(COMMANDS2)
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_1.0_v2.txt ../input/APCOM_loc_fine1_1.0.txt ./option_APCOM_SIGA_e.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.8_e: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.8.txt ../input/APCOM_loc_fine1_0.8.txt ./option_APCOM_SIGA_e.ini $(COMMANDS2)
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.8_v2.txt ../input/APCOM_loc_fine1_0.8.txt ./option_APCOM_SIGA_e.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.6_e: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.6.txt ../input/APCOM_loc_fine1_0.6.txt ./option_APCOM_SIGA_e.ini $(COMMANDS2)
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.6_v2.txt ../input/APCOM_loc_fine1_0.6.txt ./option_APCOM_SIGA_e.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.4_e: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.4.txt ../input/APCOM_loc_fine1_0.4.txt ./option_APCOM_SIGA_e.ini $(COMMANDS2)
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.4_v2.txt ../input/APCOM_loc_fine1_0.4.txt ./option_APCOM_SIGA_e.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.2_e: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.2.txt ../input/APCOM_loc_fine1_0.2.txt ./option_APCOM_SIGA_e.ini $(COMMANDS2)
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.2_v2.txt ../input/APCOM_loc_fine1_0.2.txt ./option_APCOM_SIGA_e.ini $(COMMANDS2)

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
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_1.0.txt ../input/APCOM_loc_fine1_1.0.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_1.0_v2.txt ../input/APCOM_loc_fine1_1.0.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.8_ep: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.8.txt ../input/APCOM_loc_fine1_0.8.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.8_v2.txt ../input/APCOM_loc_fine1_0.8.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.6_ep: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.6.txt ../input/APCOM_loc_fine1_0.6.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.6_v2.txt ../input/APCOM_loc_fine1_0.6.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.4_ep: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.4.txt ../input/APCOM_loc_fine1_0.4.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.4_v2.txt ../input/APCOM_loc_fine1_0.4.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.2_ep: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.2.txt ../input/APCOM_loc_fine1_0.2.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.2_v2.txt ../input/APCOM_loc_fine1_0.2.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_1.0_ep_nb: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_1.0.txt ../input/APCOM_loc_fine1_1.0.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_1.0_v2.txt ../input/APCOM_loc_fine1_1.0.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.8_ep_nb: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.8.txt ../input/APCOM_loc_fine1_0.8.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.8_v2.txt ../input/APCOM_loc_fine1_0.8.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.6_ep_nb: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.6.txt ../input/APCOM_loc_fine1_0.6.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.6_v2.txt ../input/APCOM_loc_fine1_0.6.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.4_ep_nb: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.4.txt ../input/APCOM_loc_fine1_0.4.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.4_v2.txt ../input/APCOM_loc_fine1_0.4.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.2_ep_nb: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.2.txt ../input/APCOM_loc_fine1_0.2.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.2_v2.txt ../input/APCOM_loc_fine1_0.2.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_1.0_ep_10gp: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_1.0.txt ../input/APCOM_loc_fine1_1.0.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_1.0_v2.txt ../input/APCOM_loc_fine1_1.0.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.8_ep_10gp: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.8.txt ../input/APCOM_loc_fine1_0.8.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.8_v2.txt ../input/APCOM_loc_fine1_0.8.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.6_ep_10gp: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.6.txt ../input/APCOM_loc_fine1_0.6.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.6_v2.txt ../input/APCOM_loc_fine1_0.6.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.4_ep_10gp: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.4.txt ../input/APCOM_loc_fine1_0.4.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.4_v2.txt ../input/APCOM_loc_fine1_0.4.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.2_ep_10gp: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.2.txt ../input/APCOM_loc_fine1_0.2.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.2_v2.txt ../input/APCOM_loc_fine1_0.2.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_1.0_ep_nb_10gp: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_1.0.txt ../input/APCOM_loc_fine1_1.0.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_1.0_v2.txt ../input/APCOM_loc_fine1_1.0.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.8_ep_nb_10gp: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.8.txt ../input/APCOM_loc_fine1_0.8.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.8_v2.txt ../input/APCOM_loc_fine1_0.8.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.6_ep_nb_10gp: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.6.txt ../input/APCOM_loc_fine1_0.6.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.6_v2.txt ../input/APCOM_loc_fine1_0.6.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.4_ep_nb_10gp: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.4.txt ../input/APCOM_loc_fine1_0.4.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.4_v2.txt ../input/APCOM_loc_fine1_0.4.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.2_ep_nb_10gp: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.2.txt ../input/APCOM_loc_fine1_0.2.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.2_v2.txt ../input/APCOM_loc_fine1_0.2.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_1.0_ep_linear_10gp: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_1.0.txt ../input/APCOM_loc_fine1_1.0.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_1.0_v2.txt ../input/APCOM_loc_fine1_1.0.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.8_ep_linear_10gp: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.8.txt ../input/APCOM_loc_fine1_0.8.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.8_v2.txt ../input/APCOM_loc_fine1_0.8.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.6_ep_linear_10gp: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.6.txt ../input/APCOM_loc_fine1_0.6.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.6_v2.txt ../input/APCOM_loc_fine1_0.6.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.4_ep_linear_10gp: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.4.txt ../input/APCOM_loc_fine1_0.4.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.4_v2.txt ../input/APCOM_loc_fine1_0.4.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.2_ep_linear_10gp: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.2.txt ../input/APCOM_loc_fine1_0.2.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.2_v2.txt ../input/APCOM_loc_fine1_0.2.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_1.0_ep_cos_10gp: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_1.0.txt ../input/APCOM_loc_fine1_1.0.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_1.0_v2.txt ../input/APCOM_loc_fine1_1.0.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.8_ep_cos_10gp: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.8.txt ../input/APCOM_loc_fine1_0.8.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.8_v2.txt ../input/APCOM_loc_fine1_0.8.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.6_ep_cos_10gp: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.6.txt ../input/APCOM_loc_fine1_0.6.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.6_v2.txt ../input/APCOM_loc_fine1_0.6.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.4_ep_cos_10gp: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.4.txt ../input/APCOM_loc_fine1_0.4.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.4_v2.txt ../input/APCOM_loc_fine1_0.4.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

APCOM_SIGA_fine1_0.2_ep_cos_10gp: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.2.txt ../input/APCOM_loc_fine1_0.2.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.2_v2.txt ../input/APCOM_loc_fine1_0.2.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

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

paper_hole_IGA: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_paper_hole_IGA.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_IGA.txt ./option_paper_hole_IGA.ini $(COMMANDS2)

paper_hole_IGA_200: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_paper_hole_IGA.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_IGA.txt ./option_paper_hole_IGA.ini $(COMMANDS2)

paper_hole_SIGA_conv10: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_paper_hole_SIGA.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc_oldBC.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)

paper_hole_SIGA_e_nb: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_paper_hole_SIGA.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc_oldBC.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)

paper_hole_SIGA_linear150_aux: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_paper_hole_SIGA.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc_oldBC.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)

paper_hole_SIGA_b150_aux: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_paper_hole_SIGA.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc_oldBC.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)

paper_hole_SIGA_cos150_aux: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_paper_hole_SIGA.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc_oldBC.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)

paper_hole_SIGA_nb200: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_paper_hole_SIGA.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc_oldBC.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)

paper_hole_SIGA_linear200: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_paper_hole_SIGA.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc_oldBC.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)

paper_hole_SIGA_b200: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_paper_hole_SIGA.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc_oldBC.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)

paper_hole_SIGA_cos200: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_paper_hole_SIGA.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc_oldBC.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)

paper_hole_SIGA_linear200_aux: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_paper_hole_SIGA.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc_oldBC.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)

paper_hole_SIGA_b200_aux: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_paper_hole_SIGA.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc_oldBC.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)

paper_hole_SIGA_cos200_aux: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_paper_hole_SIGA.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc_oldBC.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)
# 	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)

paper_hole_fine_IGA_200: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_paper_hole_IGA.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_IGA2.txt ./option_paper_hole_IGA.ini $(COMMANDS2)

paper_hole_fine_SIGA_nb200: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_paper_hole_SIGA2.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc3.txt ./option_paper_hole_SIGA2.ini $(COMMANDS2)

paper_hole_fine_SIGA_linear200: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_paper_hole_SIGA2.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc3.txt ./option_paper_hole_SIGA2.ini $(COMMANDS2)

paper_hole_fine_SIGA_b200: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_paper_hole_SIGA2.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc3.txt ./option_paper_hole_SIGA2.ini $(COMMANDS2)

paper_hole_fine_SIGA_cos200: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_paper_hole_SIGA2.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc3.txt ./option_paper_hole_SIGA2.ini $(COMMANDS2)

paper_hole_SIGA_vtk: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_paper_hole_SIGA.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)

debug_overlap: compile rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_debug_overlap.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src/$(PROGRAM_NAME).x ../input/debug_overlap_glo.txt ../input/debug_overlap_loc.txt ./option_debug_overlap.ini $(COMMANDS2)

debug_cylinder_B4: compile rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_debug_cylinder_B4.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src/$(PROGRAM_NAME).x ../input/debug_cylinder_B4.txt ./option_debug_cylinder_B4.ini $(COMMANDS2)

check_BC: compile rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_check_BC.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src/$(PROGRAM_NAME).x ../input/check_BC_glo.txt ../input/check_BC_loc.txt ./option_check_BC.ini $(COMMANDS2)

bf_contour: compile rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_bf_contour.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src/$(PROGRAM_NAME).x ../input/bf_contour_glo.txt ../input/bf_contour_loc.txt ./option_bf_contour.ini $(COMMANDS2)

bf_contour_v2: compile rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_bf_contour_v2.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src/$(PROGRAM_NAME).x ../input/bf_contour_glo_v2.txt ../input/bf_contour_loc_v2.txt ./option_bf_contour_v2.ini $(COMMANDS2)

# Hole problem, 200MPa, 10gp
paperHole_IGA_200: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_paper_hole_IGA.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_IGA.txt ./option_paper_hole_IGA.ini $(COMMANDS2)

paperHole_SIGA_200_consistentInt: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_paper_hole_SIGA.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc_oldBC.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)

paperHole_SIGA_200_cubic: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_paper_hole_SIGA.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc_oldBC.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)

paperHole_SIGA_200_cubic_aux: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_paper_hole_SIGA.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc_oldBC.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)

paperHole_SIGA_200_linear: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_paper_hole_SIGA.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc_oldBC.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)

paperHole_SIGA_200_linear_aux: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_paper_hole_SIGA.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc_oldBC.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)

paperHole_SIGA_200_cos: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_paper_hole_SIGA.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc_oldBC.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)

paperHole_SIGA_200_cos_aux: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_paper_hole_SIGA.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/paper_hole_SIGA_glo.txt ../input/paper_hole_SIGA_loc_oldBC.txt ./option_paper_hole_SIGA.ini $(COMMANDS2)

# IGA, 200MPa, 10gp
paperCrack_IGA_1.0_ep_200: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_IGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_iga_fine1_1.0.txt ./option_APCOM_IGA_ep.ini $(COMMANDS2)

paperCrack_IGA_0.8_ep_200: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_IGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_iga_fine1_0.8.txt ./option_APCOM_IGA_ep.ini $(COMMANDS2)

paperCrack_IGA_0.6_ep_200: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_IGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_iga_fine1_0.6.txt ./option_APCOM_IGA_ep.ini $(COMMANDS2)

paperCrack_IGA_0.4_ep_200: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_IGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_iga_fine1_0.4.txt ./option_APCOM_IGA_ep.ini $(COMMANDS2)

paperCrack_IGA_0.2_ep_200: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_IGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_iga_fine1_0.2.txt ./option_APCOM_IGA_ep.ini $(COMMANDS2)

# SIGA, consistent integration, 200MPa, 10gp
paperCrack_SIGA_1.0_ep_200_consistentInt: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_1.0_v2.txt ../input/APCOM_loc_fine1_1.0.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

paperCrack_SIGA_0.8_ep_200_consistentInt: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.8_v2.txt ../input/APCOM_loc_fine1_0.8.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

paperCrack_SIGA_0.6_ep_200_consistentInt: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.6_v2.txt ../input/APCOM_loc_fine1_0.6.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

paperCrack_SIGA_0.4_ep_200_consistentInt: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.4_v2.txt ../input/APCOM_loc_fine1_0.4.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

paperCrack_SIGA_0.2_ep_200_consistentInt: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.2_v2.txt ../input/APCOM_loc_fine1_0.2.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

# SIGA, cubic blending, 200MPa, 10gp
paperCrack_SIGA_1.0_ep_200_cubic: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_1.0_v2.txt ../input/APCOM_loc_fine1_1.0.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

paperCrack_SIGA_1.0_ep_200_cubic_5and10gp: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_1.0_v2.txt ../input/APCOM_loc_fine1_1.0.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

paperCrack_SIGA_1.0_ep_200_cubic_aux: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_1.0_v2.txt ../input/APCOM_loc_fine1_1.0.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

paperCrack_SIGA_0.8_ep_200_cubic: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.8_v2.txt ../input/APCOM_loc_fine1_0.8.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

paperCrack_SIGA_0.6_ep_200_cubic: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.6_v2.txt ../input/APCOM_loc_fine1_0.6.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

paperCrack_SIGA_0.4_ep_200_cubic: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.4_v2.txt ../input/APCOM_loc_fine1_0.4.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

paperCrack_SIGA_0.2_ep_200_cubic: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.2_v2.txt ../input/APCOM_loc_fine1_0.2.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

# SIGA, linear blending, 200MPa, 10gp
paperCrack_SIGA_1.0_ep_200_linear: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_1.0_v2.txt ../input/APCOM_loc_fine1_1.0.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

paperCrack_SIGA_1.0_ep_200_linear_aux: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_1.0_v2.txt ../input/APCOM_loc_fine1_1.0.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

paperCrack_SIGA_0.8_ep_200_linear: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.8_v2.txt ../input/APCOM_loc_fine1_0.8.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

paperCrack_SIGA_0.6_ep_200_linear: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.6_v2.txt ../input/APCOM_loc_fine1_0.6.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

paperCrack_SIGA_0.4_ep_200_linear: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.4_v2.txt ../input/APCOM_loc_fine1_0.4.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

paperCrack_SIGA_0.2_ep_200_linear: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.2_v2.txt ../input/APCOM_loc_fine1_0.2.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

# SIGA, cosine blending, 200MPa, 10gp
paperCrack_SIGA_1.0_ep_200_cos: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_1.0_v2.txt ../input/APCOM_loc_fine1_1.0.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

paperCrack_SIGA_1.0_ep_200_cos_aux: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_1.0_v2.txt ../input/APCOM_loc_fine1_1.0.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

paperCrack_SIGA_0.8_ep_200_cos: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.8_v2.txt ../input/APCOM_loc_fine1_0.8.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

paperCrack_SIGA_0.6_ep_200_cos: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.6_v2.txt ../input/APCOM_loc_fine1_0.6.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

paperCrack_SIGA_0.4_ep_200_cos: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.4_v2.txt ../input/APCOM_loc_fine1_0.4.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

paperCrack_SIGA_0.2_ep_200_cos: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_ep.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.2_v2.txt ../input/APCOM_loc_fine1_0.2.txt ./option_APCOM_SIGA_ep.ini $(COMMANDS2)

# IGA, 200MPa, 10gp
paperCrack_IGA_1.0_e_200: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_IGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_iga_fine1_1.0.txt ./option_APCOM_IGA_e.ini $(COMMANDS2)

paperCrack_IGA_0.8_e_200: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_IGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_iga_fine1_0.8.txt ./option_APCOM_IGA_e.ini $(COMMANDS2)

paperCrack_IGA_0.6_e_200: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_IGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_iga_fine1_0.6.txt ./option_APCOM_IGA_e.ini $(COMMANDS2)

paperCrack_IGA_0.4_e_200: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_IGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_iga_fine1_0.4.txt ./option_APCOM_IGA_e.ini $(COMMANDS2)

paperCrack_IGA_0.2_e_200: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_IGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_iga_fine1_0.2.txt ./option_APCOM_IGA_e.ini $(COMMANDS2)

# SIGA, consistent integration, 200MPa, 10gp
paperCrack_SIGA_1.0_e_200_consistentInt: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_1.0_v2.txt ../input/APCOM_loc_fine1_1.0.txt ./option_APCOM_SIGA_e.ini $(COMMANDS2)

paperCrack_SIGA_0.8_e_200_consistentInt: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.8_v2.txt ../input/APCOM_loc_fine1_0.8.txt ./option_APCOM_SIGA_e.ini $(COMMANDS2)

paperCrack_SIGA_0.6_e_200_consistentInt: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.6_v2.txt ../input/APCOM_loc_fine1_0.6.txt ./option_APCOM_SIGA_e.ini $(COMMANDS2)

paperCrack_SIGA_0.4_e_200_consistentInt: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.4_v2.txt ../input/APCOM_loc_fine1_0.4.txt ./option_APCOM_SIGA_e.ini $(COMMANDS2)

paperCrack_SIGA_0.2_e_200_consistentInt: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.2_v2.txt ../input/APCOM_loc_fine1_0.2.txt ./option_APCOM_SIGA_e.ini $(COMMANDS2)

# SIGA, cubic blending, 200MPa, 10gp
paperCrack_SIGA_1.0_e_200_cubic: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_1.0_v2.txt ../input/APCOM_loc_fine1_1.0.txt ./option_APCOM_SIGA_e.ini $(COMMANDS2)

paperCrack_SIGA_0.8_e_200_cubic: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.8_v2.txt ../input/APCOM_loc_fine1_0.8.txt ./option_APCOM_SIGA_e.ini $(COMMANDS2)

paperCrack_SIGA_0.6_e_200_cubic: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.6_v2.txt ../input/APCOM_loc_fine1_0.6.txt ./option_APCOM_SIGA_e.ini $(COMMANDS2)

paperCrack_SIGA_0.4_e_200_cubic: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.4_v2.txt ../input/APCOM_loc_fine1_0.4.txt ./option_APCOM_SIGA_e.ini $(COMMANDS2)

paperCrack_SIGA_0.2_e_200_cubic: compile_optimized rm mkdir_output
	cd ./analysis/material && ./cp.sh ./option_APCOM_SIGA_e.ini
	cd ./analysis/$(OUTPUT_DIR) && . ../setenv.sh && $(COMMANDS1) ../../src_fast/$(PROGRAM_NAME).x ../input/APCOM_glo_fine1_0.2_v2.txt ../input/APCOM_loc_fine1_0.2.txt ./option_APCOM_SIGA_e.ini $(COMMANDS2)

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