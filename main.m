rho=0.001;

belief_samples=getbelief(rho);


[hyper_f_shed,hyper_f_unshed ]= pesues(rand(3,14)*30,rand(3,14)*30,belief_samples);


while 1
    [hyper_f_shed,hyper_f_unshed ]= pesues(hyper_f_shed,hyper_f_unshed ,belief_samples)
end
