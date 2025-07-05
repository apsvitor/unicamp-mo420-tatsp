# Define o diretório das instâncias
# INSTANCES_DIR = ../../instances/instances_release_2
INSTANCES_DIR = instances

# Encontra todos os arquivos .txt no diretório de instâncias
# A função $(wildcard ...) expande para uma lista de arquivos que correspondem ao padrão
# A função $(notdir ...) remove o caminho do diretório, deixando apenas o nome do arquivo
# A função $(basename ...) remove a extensão do arquivo
INSTANCES := $(basename $(notdir $(wildcard $(INSTANCES_DIR)/*.txt)))

TEST_INSTANCE := $(INSTANCES_DIR)/inst-slide8.txt

.PHONY: all test clean

# Regra 'all': executa o problema para cada instância
all: $(foreach instance, $(INSTANCES), $(instance).log)

# Regra de pattern matching para gerar cada arquivo de log
# $@ é o nome do target (ex: grf101.log)
# $* é o nome do target sem a extensão (ex: grf101)
%.log:
	@echo "Executando para a instância: $*"
	julia --project=. trigger_arc_tsp.jl \
		--inputfilename $(INSTANCES_DIR)/$*.txt \
		--seednumber 1234 \
		--maxtime_lb_lp 1010 \
		--maxtime_lb_rlxlag 1020 \
		--maxtime_lb_colgen 1030 \
		--maxtime_ub_lp 1040 \
		--maxtime_ub_rlxlag 1050 \
		--maxtime_ub_colgen 1060 \
		--maxtime_ilp 1070 \
		--ra 291207 \
		--logfile logs/test_$*.log \
		--verbose

test:
	@echo Execução de instância de teste:
	julia --project=. trigger_arc_tsp.jl \
			--inputfilename $(TEST_INSTANCE) \
			--seednumber 1234 \
			--maxtime_lb_lp 1010 \
			--maxtime_lb_rlxlag 1020 \
			--maxtime_lb_colgen 1030 \
			--maxtime_ub_lp 1040 \
			--maxtime_ub_rlxlag 1050 \
			--maxtime_ub_colgen 1060 \
			--maxtime_ilp 1070 \
			--ra 291207 \
			--logfile trigger_arc_tsp.log \
			--verbose

# Regra 'clean': remove todos os arquivos de log gerados
clean:
	@echo "Removendo arquivos de log..."
	rm -f $(INSTANCES)/logs/*.log