 counter = 0 ;


        for m =1:1:sample_number

            %calculate J_k for all samples
            J_k= dot(belief_samples(:,m),hyper_planes(1,:));
            if no_planes >1
                for k=2:1:no_planes
                    if dot(belief_samples(:,m),hyper_planes(k,:)) <= J_k
                        J_k = dot(belief_samples(:,m),hyper_planes(k,:));
                    end

                end
            end

            %calculate J_k+1 for all samples

            J_k1= dot(belief_samples(:,m),hyper_planes_future(1,:));
            if no_planes_future>1
                for k=2:1:no_planes_future
                    if dot(belief_samples(:,m),hyper_planes_future(k,:)) <= J_k1
                        J_k1 = dot(belief_samples(:,m),hyper_planes_future(k,:));
                    end

                end
            end

            if J_k1>J_k
                counter = counter+1;
            end

        end

        if counter ==0 && size_hyper_planes_non_shed_future(2)~=1 &&size_hyper_planes_shed_future(2)~=1
                break
        end