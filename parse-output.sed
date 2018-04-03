/^[[:space:]]*[[:digit:]]\+$/ {
N
s_\n[^-]*-\+_/_
N
s/\n//
s/[[:space:]]*//g
p
}
