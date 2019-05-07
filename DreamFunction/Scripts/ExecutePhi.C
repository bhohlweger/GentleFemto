#include "GetCorrelationsPhi.C"

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";
  GetCorrelationsPhi(filename, "Results", prefix, addon);
  return 1;
}
