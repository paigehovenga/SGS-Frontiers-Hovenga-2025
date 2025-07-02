function [pdf_out, CI05, CI95, CI025, CI975] = estimate_pdf_CI(data, xx);

% Description:
% Computes the confidence intervals from SGS synthetic data or other dataset

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inputs:
% data = matrix of the bootstrapped SGS synthetic time series
% xx = values to evalues the PDF at

% Outputs:
% pdf_out = sorted PDFs
% CI05, CI95, CI025, CI975 = confidence intervals 

% Provenance: 
% Affiliation: Physical Sciences Laboratory, NOAA and Cooperative Institue for Research in
% Environmental Sciences, CU Boulder
% Written: ~Feb 2017 by Matt Newman
% Modified: Jun 2025 by Paige Hovenga (contact paige.hovenga@noaa.gov)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pdf_out = [];

    %Create PDF using ksdensity
    for i = 1:size(data,1);
        [tmp1,x]= ksdensity(data(i,:), xx); %could create smoother
        pdf_out = cat(1,pdf_out, tmp1);
    end

    %Confidence Intervals
    for i=1:size(pdf_out,2);
        pdfsort(i,:)=sort(pdf_out(:,i));
    end

    pdf_out = pdfsort;

    if size(pdf_out,2) > 1; 
        CI05=pdf_out(:,(round(size(pdf_out,2)*0.05)));
        CI95=pdf_out(:,size(pdf_out,2)-round((size(pdf_out,2)*0.05)-1));
        CI025=pdf_out(:,round(size(pdf_out,2)*0.025));
        CI975=pdf_out(:,round(size(pdf_out,2)-(size(pdf_out,2)*0.025-1)));
    else %if bootstapping done on outside of code
        CI05 = NaN;
        CI95 = NaN;
        CI025 = NaN;
        CI975 = NaN; 
    end

end