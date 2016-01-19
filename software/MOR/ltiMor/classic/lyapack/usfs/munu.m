function usfs = munu
  usfs.pre       = @munu_pre;
  usfs.post      = @munu_pst;

  usfs.mi        = @munu_m_i;
  usfs.m         = @munu_m;
  usfs.md        = @munu_m_d;

  usfs.li        = @munu_l_i;
  usfs.l         = @munu_l;
  usfs.ld        = @munu_l_d;
  
  usfs.si        = @munu_s_i;
  usfs.s         = @munu_s;
  usfs.sd        = @munu_s_d;
  
  usfs.e         = @munu_e;
  