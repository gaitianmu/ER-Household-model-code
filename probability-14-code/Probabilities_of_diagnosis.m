[A,B,C]=xlsread([cd '/Covid19CasesWH0706.xlsx']);
Onset_data=C;
Onset_data(1,:)=[];
num_days=size(Onset_data);
num_days=num_days(1);

Onset_vector=zeros(0,num_days);
for day=1:num_days
    Onset_vector(day)=Onset_data{day,2};
end


[A,B,C]=xlsread([cd '/Phasesro.xlsx']);
Phase_data=C;
Phase_data(1,:)=[];

num_phases=size(Phase_data);
num_phases=num_phases(1);

Dq_vector=zeros(1,num_phases);
for i=1:num_phases
    Dq_vector(i)=Phase_data{i,2};
end


changing_points_vector=zeros(1,num_phases);

for i=1:num_phases
    date_target=Phase_data{i,1};
    for j=1:num_days
        date_temp=Onset_data{j,1};
        if length(date_target)==length(date_temp)&& prod(date_target==date_temp)
            changing_points_vector(i)=j;
            break;
        end
    end
end

D_i_vector=[2.9 5 7 8 9 10 11 12 13 14];
num_cases=length(D_i_vector);

for case_temp=1:num_cases
    Diagnosis_probability_vector=zeros(1,num_days);
    D_i=D_i_vector(case_temp); 
    num_samples=10^5; %number of samples in approximating probability
    
    
    for day=1:num_days
        [D_i_vector(case_temp) day]
        num_diagnosis=0;
        for sample=1:num_samples
            day_recover=day+exprnd(D_i);
            day_diagnosis=Diagnosis_date(day,changing_points_vector,Dq_vector);
            if day_diagnosis<day_recover
                num_diagnosis=num_diagnosis+1;
            end
        end
        p_diagnosis=num_diagnosis/num_samples;
        Diagnosis_probability_vector(day)=p_diagnosis;
    end
    
    date_label_vector=cell(1,num_days);
    date_start=Onset_data{1,1};
    
    for day=1:num_days
        day_temp = datetime(datestr(date_start,'yyyy-mm-dd'))+caldays(day-1);  %%%%-2
        day_temp=datestr(day_temp,'mm-dd');
        date_label_vector{day}=day_temp;
    end
    
    
    Supposed_diagnosis_data=cell(num_days,2);
    for i=1:num_days
        Supposed_diagnosis_data(i,1)=date_label_vector(i);
        Supposed_diagnosis_data(i,2)={Diagnosis_probability_vector(i)};
    end
    %xlswrite([cd '/Appoximated_diagnosis_probability_' num2str(D_i) '.xlsx'],Supposed_diagnosis_data);
    writecell(Supposed_diagnosis_data,['ro_0706_Appoximated_diagnosis_probability_ ' num2str(D_i) '.xlsx'])
    
end

%plot(Diagnosis_probability_vector);
%set(gca, 'XTick', 0:5:num_days-1, 'XTickLabel', date_label_vector(1:5:num_days))


