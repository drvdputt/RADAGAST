set_property(TARGET ${TARGET} PROPERTY CXX_STANDARD 14)

target_compile_options(${TARGET} PRIVATE -Wshadow -Wall -Wextra -Wno-missing-braces -Wno-sign-compare -Werror=return-type -pedantic)

