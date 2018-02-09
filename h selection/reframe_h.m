ha_100=hvec1(index_a1);
ha_200=hvec1(index_a2);
ha_300=hvec2(index_a3);
ha_400=hvec2(index_a4);

hb_100=hvec(index_b1);
hb_200=hvec(index_b2);
hb_300=hvec(index_b3);
hb_400=hvec(index_b4);

hc_100=hvec1(index_c1);
hc_200=hvec2(index_c2);
hc_300=hvec2(index_c3);
hc_400=hvec2(index_c4);

save h100 ha_100 hb_100 hc_100 hrho_100
save h200 ha_200 hb_200 hc_200 hrho_200
save h300 ha_300 hb_300 hc_300 hrho_300
save h400 ha_400 hb_400 hc_400 hrho_400

load h_sim300_mn100new
cv=cvb; %n*24 plot the lscv vs h curve
plot(cv'-kron(ones(24,1),mean(cv')))


[ise_a,index_a]=min(cva');
[ise_b,index_b]=min(cvb');
[ise_c,index_c]=min(cvc');
[ise_rho,index_rho]=min(cvrho');
hvec=0.0025:0.0025:0.06;
ha_200=hvec(index_a);
hb_200=hvec(index_b);
hc_200=hvec(index_c);
hrho_200=hvec(index_rho);
mean([ha_200' hb_200' hc_200' hrho_200'])
save h200 ha_200 hb_200 hc_200 hrho_200

load h100
mean([ha_400' hb_400' hc_400' hrho_400'])
boxplot([ha_300' hb_300' hc_300' hrho_300'])
boxplot([ha_200' hb_200' hc_200' hrho_200'])
mean([ha_100' hb_100' hc_100' hrho_100'])





