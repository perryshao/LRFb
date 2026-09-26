/* Validation-only scalar build: arm64 has no x86 CPU feature flags. */
static inline void _vl_cpuid(int *info, int function)
{
    (void)function;
    info[0] = info[1] = info[2] = info[3] = 0;
}
