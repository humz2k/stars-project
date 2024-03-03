RAYLIB_DIR ?= ../ext/raylib/src
RAYLIB_FLAGS ?= -framework CoreVideo -framework IOKit -framework Cocoa -framework GLUT -framework OpenGL

BUILD_DIR ?= build
SOURCE_DIR ?= src
INCLUDE_DIR ?= include

SOURCES := $(shell find $(SOURCE_DIR) -name '*.c')
OBJECTS := $(SOURCES:%.c=%.o)
OUTPUTS := $(OBJECTS:src%=build%) $(SCRIPT_OUTPUTS)

HEADERS := $(shell find $(SOURCE_DIR) -name '*.h')

main: driver

.PHONY: $(RAYLIB_DIR)/libraylib.a
$(RAYLIB_DIR)/libraylib.a:
	cd $(RAYLIB_DIR) && $(MAKE) MACOSX_DEPLOYMENT_TARGET=10.9

driver: $(OUTPUTS) $(RAYLIB_DIR)/libraylib.a
	gcc-13 -O3 $^ -o $@ $(RAYLIB_FLAGS) -fopenmp

build/%.o: $(SOURCE_DIR)/%.c $(HEADERS) | $(BUILD_DIR)
	gcc-13 -O3 -I$(RAYLIB_DIR) -I$(SOURCE_DIR) -I$(INCLUDE_DIR) -c -o $@ $< -fopenmp

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf build
	rm -rf driver