#define sha2_starts oracle_internal_sha2_starts
#define sha2_update oracle_internal_sha2_update
#define sha2_finish oracle_internal_sha2_finish
#define sha2_file oracle_internal_sha2_file
#define sha2_csum oracle_internal_sha2_csum
#define sha2_hmac oracle_internal_sha2_hmac
#define sha2_self_test oracle_internal_sha2_self_test

#include "sha2.c"

void oracle_sha2_starts(sha2_context *ctx)
{
    oracle_internal_sha2_starts(ctx);
}

void oracle_sha2_process(sha2_context *ctx, unsigned char data[64])
{
    sha2_process(ctx, data);
}

void oracle_sha2_update(sha2_context *ctx, unsigned char *input, int ilen)
{
    oracle_internal_sha2_update(ctx, input, ilen);
}

void oracle_sha2_finish(sha2_context *ctx, unsigned char output[32])
{
    oracle_internal_sha2_finish(ctx, output);
}

void oracle_sha2_csum(unsigned char *input, int ilen, unsigned char output[32])
{
    oracle_internal_sha2_csum(input, ilen, output);
}
