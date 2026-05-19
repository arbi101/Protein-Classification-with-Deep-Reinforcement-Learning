from tests.test_hc import generate_2d_structure as run_hc
from tests.test_sa import generate_2d_structure_sa as run_sa
from tests.test_remc import sequence_to_hp

p2 = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
hp_str = sequence_to_hp(p2)

_, en_hc = run_hc(hp_str, iterations=300000)
_, en_sa = run_sa(hp_str, iterations=300000, initial_t=30.0, final_t=0.001)

print(f"HC: {en_hc}")
print(f"SA: {en_sa}")
