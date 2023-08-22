void SAIGEClass::scoreTestFast_noadjCov_multi(arma::mat & t_GMat_centered,
                     arma::vec & t_Beta,
                     arma::vec & t_seBeta,
                     std::vector <std::string> & t_pval_str,
                     arma::vec &t_Tstat,
                     arma::vec &t_var1,
                     arma::vec &t_var2,
                     arma::uvec & t_skipSPAind_vec,
                     arma::vec & t_pval_noadj_vec)
{
t_Tstat = t_GMat_centered.t() * m_res_sample / (m_tauvec[0]);
//arma::mat mu2Mat = arma::diag(m_mu2_sample);
 std::cout << "t_GMat here2c " << std::endl;
arma::mat t_GMat_centered_weight = (t_GMat_centered.each_col()) % m_mu2_sample;
 std::cout << "t_GMat here2d " << std::endl;
t_GMat_centered_weight = t_GMat_centered_weight % t_GMat_centered;
 std::cout << "t_GMat here2e " << std::endl;
t_var2 = (arma::sum(t_GMat_centered_weight, 0).t());
std::cout << "t_GMat here2f " << std::endl;
t_var1 = t_var2 * m_varRatioVal;
 std::cout << "t_GMat here2g " << std::endl;
arma::vec stat_vec = arma::pow(t_Tstat, 2) / t_var1;
double pval, stat;
std::cout << "t_GMat here2h " << std::endl;

int numMarkers = t_var1.n_elem;

std::cout << "numMarkers " << numMarkers << std::endl;
t_pval_noadj_vec.set_size(numMarkers);

for(int i = 0; i < numMarkers; i++){
        if (t_var1(i) <= std::numeric_limits<double>::min()){
                pval = 1;
                stat = 0;
        }else{
                stat = stat_vec(i);
                pval = Rf_pchisq(stat, 1.0, 0, 0);
        }
std::cout << "t_GMat here2h0 " << std::endl;
        char pValueBuf[100];
        if (pval != 0){
                sprintf(pValueBuf, "%.6E", pval);
        }else {
                double log10p = log10(2.0) - M_LOG10E*stat/2 - 0.5*log10(stat*2*M_PI);
                int exponent = floor(log10p);
                double fraction = pow(10.0, log10p - exponent);
                if (fraction >= 9.95) {
                        fraction = 1;
                         exponent++;
                }
                sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
        }
std::cout << "t_GMat here2h1 " << std::endl;
        std::string buffAsStdStr = pValueBuf;
        t_pval_str.push_back(buffAsStdStr);
        t_pval_noadj_vec(i) = pval;
}
std::cout << "t_GMat here2f " << std::endl;
t_Beta = t_Tstat/t_var1;
t_seBeta = arma::abs(t_Beta) / sqrt(arma::abs(stat_vec));

arma::vec StdStat_vec = arma::abs(t_Tstat) / arma::sqrt(t_var1);
//t_skipSPAind_vec.elem( find(StdStat_vec > m_SPA_Cutoff) ).zeros();

t_skipSPAind_vec = arma::find(StdStat_vec > m_SPA_Cutoff);

}
