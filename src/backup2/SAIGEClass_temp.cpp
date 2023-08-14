void SAIGEClass::scoreTestFast_noadjCov_multi(arma::mat & t_GMat_centered,
                     arma::vec & t_Beta,
                     arma::vec & t_seBeta,
                     std::vector <std::string> & t_pval_str,
                     arma::vec &t_Tstat,
                     arma::vec &t_var1,
                     arma::vec &t_var2, 
		     arma::uvec & t_skipSPAind_vec, 
		     arma::vec & t_pval_noadj_vec){

t_Tstat = t_GMat_centered.t() * m_res_sample / (m_tauvec[0]);
t_var2 = t_GMat_centered.t() * arma::diag(m_mu2_sample) * t_GMat_centered;
t_var1 = t_var2 % t_varRatio_vec;
arma::vec stat_vec = arma::pow(t_Tstat, 2) / t_var1;
double pval, stat;


int numMarkers = t_var1.n_elem;

for(int i = 0; i < numMarkers; i++){
	if (t_var1(i) <= std::numeric_limits<double>::min()){
		pval = 1;
        	stat = 0;	
	}else{
		pval = Rf_pchisq(stat, 1.0, 0, 0);
		stat = stat_vec(i);
	}
	char pValueBuf[100];
        if (pval != 0)
        	sprintf(pValueBuf, "%.6E", pval);
    	else {
        	double log10p = log10(2.0) - M_LOG10E*stat/2 - 0.5*log10(stat*2*M_PI);
        	int exponent = floor(log10p);
        	double fraction = pow(10.0, log10p - exponent);
        	if (fraction >= 9.95) {
          		fraction = 1;
          		 exponent++;
         	}
        	sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
    	}
	std::string buffAsStdStr = pValueBuf;
	t_pval_str.at(i) = buffAsStdStr;
	t_pval_noadj_vec(i) = pval;
}
t_Beta = t_Tstat/t_var1;
t_seBeta = fabs(t_Beta) / sqrt(fabs(stat_vec));

arma::vec StdStat_vec = arma::abs(t_Tstat) / arma::sqrt(t_var1);
//t_skipSPAind_vec.elem( find(StdStat_vec > m_SPA_Cutoff) ).zeros();

t_skipSPAind_vec = arma::find(StdStat_vec > m_SPA_Cutoff)

}

