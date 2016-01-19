function usfs = au
  usfs.pre       = @au_pre;
  usfs.post      = @au_pst;

  usfs.mi        = @au_m_i;
  usfs.m         = @au_m;
  usfs.md        = @au_m_d;

  usfs.li        = @au_l_i;
  usfs.l         = @au_l;
  usfs.ld        = @au_l_d;
  
  usfs.si        = @au_s_i;
  usfs.s         = @au_s;
  usfs.sd        = @au_s_d;
  
  usfs.e         = @au_e;
  
