#include <uuid/uuid.h>
#include "uuid_str.h"
using namespace std;

std::string uuid_str(void)
{
  unsigned char uuid_c[16];
  char s[36];

  // generate uuid
  uuid_generate_time(uuid_c);

  // convert to string
  uuid_unparse(uuid_c,s);

  return string(s);
}
