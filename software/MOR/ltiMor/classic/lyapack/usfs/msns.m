function usfs = msns
  usfs.pre       = @msns_pre;
  usfs.post      = @msns_pst;

  usfs.mi        = @msns_m_i;
  usfs.m         = @msns_m;
  usfs.md        = @msns_m_d;

  usfs.li        = @msns_l_i;
  usfs.l         = @msns_l;
  usfs.ld        = @msns_l_d;
  
  usfs.si        = @msns_s_i;
  usfs.s         = @msns_s;
  usfs.sd        = @msns_s_d;
  
  usfs.e         = @msns_e;
  