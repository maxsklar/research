a = [3 1 2];

for n = 1:2;
  M = (2^n);
  data = polya_sample(a, repmat(M,1,1000));
  
  tic;
  polya_fit(data);
  time1 = toc;
  
  tic;
  polya_fit_simple(data)
  time2 = toc;

  tic;
  polya_fit_modified(data);
  time3 = toc;

  tic;
  polya_fit_simple_modified(data)
  time3 = toc;
  
  %disp([num2str(M) 9 num2str(time1) 9 num2str(time2) 9 num2str(time3)]);
end;