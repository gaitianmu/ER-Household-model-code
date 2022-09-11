function day_diagnosis=Diagnosis_date(day_onset,changing_points_vector,Dq_vector)
phase_current=find(changing_points_vector>day_onset);
if length(phase_current)==0 % if already in the last phase
    change_point_last=length(changing_points_vector);
    Dq_last=length(Dq_vector);
    day_diagnosis_temp=day_onset+exprnd(Dq_vector(Dq_last));
else
    phase_current=phase_current(1)-1;
    Dq_vector_current=Dq_vector;
    Dq_vector_current(1:phase_current-1)=[];
    changing_points_vector_current=changing_points_vector;
    changing_points_vector_current(1:phase_current)=[];
    num_changes=length(changing_points_vector_current);
    
    step=1;
    changing_point_temp=changing_points_vector_current(1);
    day_diagnosis_temp=day_onset;
    
    
    while 1
        day_diagnosis_temp=day_diagnosis_temp+exprnd(Dq_vector_current(step));
        if step>num_changes||day_diagnosis_temp<=changing_point_temp % if now in the last phase or confirmed before the end of this phase
            break;
        else %if not confirmed until the end of this phase move to the next and change diagnosis rate 
            day_diagnosis_temp=changing_point_temp;
            step=step+1;
            if step<=num_changes
                changing_point_temp=changing_points_vector_current(step);
            end
        end
    end
end
day_diagnosis=day_diagnosis_temp;
end


