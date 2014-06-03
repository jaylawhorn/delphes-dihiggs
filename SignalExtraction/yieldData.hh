#ifndef YIELD_DATA_HH
#define YIELD_DATA_HH

struct yieldData
{
  Float_t total_yield;

  Float_t total_tt, total_mt, total_et, total_em;

  Float_t tt_2, mt_2, et_2, em_2;
  Float_t tt_1, mt_1, et_1, em_1;
  Float_t tt_0, mt_0, et_0, em_0;

  Float_t tt_2_boo, mt_2_boo, et_2_boo, em_2_boo;
  Float_t tt_1_boo, mt_1_boo, et_1_boo, em_1_boo;
  Float_t tt_0_boo, mt_0_boo, et_0_boo, em_0_boo;

  Float_t tt_2_nb, mt_2_nb, et_2_nb, em_2_nb;
  Float_t tt_1_nb, mt_1_nb, et_1_nb, em_1_nb;
  Float_t tt_0_nb, mt_0_nb, et_0_nb, em_0_nb;


};

#endif
