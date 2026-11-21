# ============================================================
# Makefile для проекта "projectile" (GSL + gnuplot)
# Совместим с GNU Make >= 4.2.1
# ============================================================

# --- Идентификация проекта ----------------------------------
TARGET   := projectile
SRC      := projectile.c
OBJ      := $(SRC:.c=.o)
DEP      := $(SRC:.c=.d)
DATAFILE := trajectory.dat

# --- Инструментарий -----------------------------------------
CC       := gcc
PKGCFG   := pkg-config

# --- Флаги компиляции и линковки ----------------------------
# -MMD -MP: генерация .d-файлов с зависимостями от заголовков
# -Wall -Wextra -O2: предупреждения + оптимизация
CFLAGS   := -Wall -Wextra -O2 -std=c11 -MMD -MP
CFLAGS   += $(shell $(PKGCFG) --cflags gsl)

LDFLAGS  := $(shell $(PKGCFG) --libs-only-L gsl)
LDLIBS   := $(shell $(PKGCFG) --libs-only-l gsl) -lm

# --- Цели по умолчанию --------------------------------------
.DEFAULT_GOAL := all

.PHONY: all run plot clean distclean help

all: $(TARGET)

# --- Линковка -----------------------------------------------
$(TARGET): $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# --- Компиляция (.c -> .o) ----------------------------------
# Правило уже описано неявно, но явно задаём зависимости от заголовков
%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

# --- Подтягивание автоматически сгенерированных зависимостей
-include $(DEP)

# --- Запуск программы с сохранением данных ------------------
run: $(TARGET)
	./$(TARGET) > $(DATAFILE)

# --- Визуализация траектории через gnuplot ------------------
plot: run
	@command -v gnuplot >/dev/null 2>&1 || { \
		echo "Ошибка: gnuplot не найден. Установите: sudo apt install gnuplot"; \
		exit 1; \
	}
	gnuplot -p -e " \
		set xlabel 'Дальность, м'; \
		set ylabel 'Высота, м'; \
		set title 'Траектория с учётом сопротивления воздуха'; \
		set grid; \
		plot '$(DATAFILE)' using 2:3 with lines lw 2 title 'Траектория'"

# --- Очистка ------------------------------------------------
clean:
	$(RM) $(OBJ) $(DEP) $(TARGET)

distclean: clean
	$(RM) $(DATAFILE)

# --- Справка ------------------------------------------------
help:
	@echo "Доступные цели:"
	@echo "  all       - собрать исполняемый файл (по умолчанию)"
	@echo "  run       - запустить программу, сохранить данные в $(DATAFILE)"
	@echo "  plot      - построить график траектории через gnuplot"
	@echo "  clean     - удалить объектные файлы и бинарник"
	@echo "  distclean - удалить всё, включая файлы данных"
	@echo "  help      - показать это сообщение"
