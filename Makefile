SRC = $(wildcard src/*.c)
OBJ = $(SRC:.c=.o)

LIB = libexternalmodel.a

$(LIB): $(OBJ)
	$(AR) rcs $@ $^

.PHONY: install
install: $(LIB)
	mkdir -p $(DESTDIR)$(PREFIX)/lib
	cp $(lib) $(DESTDIR)$(PREFIX)/lib/$(lib)

.PHONY: uninstall
uninstall:
	rm -f $(DESTDIR)$(PREFIX)/lib/$(lib)

.PHONY: clean
clean:
	rm -f $(OBJ) $(LIB)
