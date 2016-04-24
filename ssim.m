tic();
input_file = fopen('xaaa');
number_of_lines = fskipl(input_file, Inf);
frewind(input_file);
cells = cell(number_of_lines, 1);
for i = 1:2*number_of_lines
    cells{i} = fscanf(input_file, '%s', 1);
end

pkg load image
indexuri = cell(number_of_lines, 1);
j=0

for i = 1:2:2*number_of_lines

    if (i==1) 
    img1=imread(strcat('/home/john/Pictures/mountain/', cells{i}));
    img2=imread(strcat('/home/john/Pictures/mountain/', cells{i+1}));
    nume_fisier_anterior=cells{i};
    [mssim, ssim_map, img1_extern, img1_squared_extern, mu11_extern, mu1_sq_extern, mu1111_extern] = ssim1_pyssim_padarray(img1, img2);
    end
    if (!strncmp(cells{i}, nume_fisier_anterior, length(cells{i})) && i>1)
    img1=imread(strcat('/home/john/Pictures/mountain/', cells{i}));
    img2=imread(strcat('/home/john/Pictures/mountain/', cells{i+1}));
    [mssim, ssim_map, img1_extern, img1_squared_extern, mu11_extern, mu1_sq_extern, mu1111_extern] = ssim1_pyssim_padarray(img1, img2);
    else
    img2=imread(strcat('/home/john/Pictures/mountain/', cells{i+1}));
    [mssim, ssim_map] = ssim1_pyssim_padarray_extern(img1_extern, img2, img1_squared_extern, mu11_extern, mu1_sq_extern, mu1111_extern);
    end
    j=j+1;
    indexuri{j}=mssim;

end

toc();


