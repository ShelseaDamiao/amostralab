(function(){
/* ---------- numeric core ---------- */
function normInv(p){
  const a=[-3.969683028665376e+01,2.209460984245205e+02,-2.759285104469687e+02,1.383577518672690e+02,-3.066479806614716e+01,2.506628277459239e+00];
  const b=[-5.447609879822406e+01,1.615858368580409e+02,-1.556989798598866e+02,6.680131188771972e+01,-1.328068155288572e+01];
  const c=[-7.784894002430293e-03,-3.223964580411365e-01,-2.400758277161838e+00,-2.549732539343734e+00,4.374664141464968e+00,2.938163982698783e+00];
  const d=[7.784695709041462e-03,3.224671290700398e-01,2.445134137142996e+00,3.754408661907416e+00];
  const pl=0.02425, ph=1-pl; let q,r;
  if(p<pl){q=Math.sqrt(-2*Math.log(p));return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);}
  if(p<=ph){q=p-0.5;r=q*q;return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q/(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);}
  q=Math.sqrt(-2*Math.log(1-p));return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
}
function tInv(p,df){ // Cornish–Fisher expansion
  const z=normInv(p), z3=z*z*z, z5=z3*z*z, z7=z5*z*z;
  return z+(z3+z)/(4*df)+(5*z5+16*z3+3*z)/(96*df*df)+(3*z7+19*z5+17*z3-15*z)/(384*df*df*df);
}
const fmt=(x,dec=0)=>Number(x).toLocaleString('pt-PT',{minimumFractionDigits:dec,maximumFractionDigits:dec});
const fmtEn=(x,dec=0)=>Number(x).toLocaleString('en-GB',{minimumFractionDigits:dec,maximumFractionDigits:dec});
const ceil=Math.ceil;
function range(a,b,s){const r=[];for(let x=a;x<=b+1e-9;x+=s)r.push(+x.toFixed(4));return r;}
const pct=x=>fmt(x*100);

const CONF=[[0.90,'90%'],[0.95,'95%'],[0.99,'99%']];
const ALPHA=[[0.01,'0,01'],[0.05,'0,05'],[0.10,'0,10']];
const SIDES=[[2,'Bilateral'],[1,'Unilateral']];
const POWER=[[0.8,'80%'],[0.85,'85%'],[0.9,'90%']];
const zA=(a,s)=>normInv(1-a/s), zB=p=>normInv(p);

/* ---------- references ---------- */
const REFS={
  cochran:{c:'Cochran, 1977',f:'Cochran, W. G. (1977). <i>Sampling Techniques</i> (3.ª ed.). Wiley.',u:'Amostra para proporções e médias; correcção para população finita.',m:'Prevalência, Média'},
  kish:{c:'Kish, 1965',f:'Kish, L. (1965). <i>Survey Sampling</i>. Wiley.',u:'Efeito de desenho (deff) em amostras por conglomerados.',m:'Prevalência, Média'},
  lohr:{c:'Lohr, 2021',f:'Lohr, S. L. (2021). <i>Sampling: Design and Analysis</i> (3.ª ed.). CRC Press.',u:'Amostragem complexa, estratificação e não resposta.',m:'Prevalência, Média'},
  fleiss:{c:'Fleiss et al., 2003',f:'Fleiss, J. L., Levin, B. e Paik, M. C. (2003). <i>Statistical Methods for Rates and Proportions</i> (3.ª ed.). Wiley.',u:'Comparação de duas proporções, com grupos de tamanho diferente.',m:'Duas proporções'},
  chow:{c:'Chow et al., 2008',f:'Chow, S.-C., Shao, J. e Wang, H. (2008). <i>Sample Size Calculations in Clinical Research</i> (2.ª ed.). Chapman &amp; Hall/CRC.',u:'Comparação de duas médias; poder e erro tipo I.',m:'Duas médias'},
  cohen:{c:'Cohen, 1988',f:'Cohen, J. (1988). <i>Statistical Power Analysis for the Behavioral Sciences</i> (2.ª ed.). Erlbaum.',u:'Tamanho de efeito padronizado (d de Cohen).',m:'Duas médias, Conglomerados'},
  borm:{c:'Borm et al., 2007',f:'Borm, G. F., Fransen, J. e Lemmens, W. A. (2007). A simple sample size formula for analysis of covariance in randomized clinical trials. <i>Journal of Clinical Epidemiology</i>, 60(12), 1234–1238.',u:'Redução da amostra com ajuste pela linha de base (ANCOVA).',m:'Duas médias'},
  hayes:{c:'Hayes e Moulton, 2017',f:'Hayes, R. J. e Moulton, L. H. (2017). <i>Cluster Randomised Trials</i> (2.ª ed.). CRC Press.',u:'Ensaios aleatorizados por conglomerados; ICC.',m:'Conglomerados, MDES'},
  donner:{c:'Donner e Klar, 2000',f:'Donner, A. e Klar, N. (2000). <i>Design and Analysis of Cluster Randomization Trials in Health Research</i>. Arnold.',u:'Efeito de desenho 1 + (m − 1)·ICC.',m:'Conglomerados'},
  bloom95:{c:'Bloom, 1995',f:'Bloom, H. S. (1995). Minimum detectable effects: A simple way to report the statistical power of experimental designs. <i>Evaluation Review</i>, 19(5), 547–556.',u:'Efeito mínimo detectável (MDES) e multiplicador t.',m:'MDES'},
  bloom05:{c:'Bloom, 2005',f:'Bloom, H. S. (2005). Randomizing groups to evaluate place-based programs. Em <i>Learning More from Social Experiments</i>. Russell Sage Foundation.',u:'MDES em desenhos por conglomerados com covariáveis.',m:'Conglomerados, MDES'},
  dong:{c:'Dong e Maynard, 2013',f:'Dong, N. e Maynard, R. (2013). PowerUp!: A tool for calculating minimum detectable effect sizes and minimum required sample sizes. <i>Journal of Research on Educational Effectiveness</i>, 6(1), 24–67.',u:'Fórmulas de MDES e amostra mínima para avaliações em educação.',m:'Conglomerados, MDES'},
  duflo:{c:'Duflo et al., 2007',f:'Duflo, E., Glennerster, R. e Kremer, M. (2007). Using randomization in development economics research: A toolkit. <i>Handbook of Development Economics</i>, vol. 4, 3895–3962.',u:'Poder estatístico em avaliações de impacto em países em desenvolvimento.',m:'Conglomerados, MDES'},
  kraft:{c:'Kraft, 2020',f:'Kraft, M. A. (2020). Interpreting effect sizes of education interventions. <i>Educational Researcher</i>, 49(4), 241–253.',u:'Referências para interpretar efeitos em educação (pequeno &lt; 0,05; médio 0,05–0,20; grande ≥ 0,20 DP).',m:'MDES'}
};

/* ---------- modules ---------- */
const MODS=[
{ id:'prevalencia', k:'Inquérito', n:'Prevalência', d:'Estimar uma percentagem com margem de erro',
  lead:'Para inquéritos que querem estimar uma proporção, como a percentagem de raparigas que concluíram a 7.ª classe.',
  fields:[
    {k:'p',l:'Prevalência esperada',u:'%',def:50,min:1,max:99,step:1,range:true,h:'Sem informação prévia, use 50%: dá a amostra mais conservadora.'},
    {k:'d',l:'Margem de erro (±)',u:'p.p.',def:5,min:0.5,max:20,step:0.5,range:true,h:'Em pontos percentuais. ±5 p.p. é o padrão em inquéritos nacionais.'},
    {k:'conf',l:'Nível de confiança',sel:CONF,def:0.95},
    {k:'N',l:'Tamanho da população',u:'pessoas',def:0,min:0,step:100,h:'Deixe 0 se a população for muito maior do que a amostra.'},
    {k:'deff',l:'Efeito de desenho',u:'deff',def:2,min:1,max:6,step:0.1,range:true,h:'1 = amostragem aleatória simples. Inquéritos domiciliares por conglomerados ficam entre 1,5 e 2,5.'},
    {k:'resp',l:'Taxa de resposta esperada',u:'%',def:90,min:30,max:100,step:1,range:true}
  ],
  presets:[
    {n:'Inquérito domiciliar por conglomerados',v:{p:50,d:5,conf:0.95,N:0,deff:2,resp:90}},
    {n:'Utentes de uma unidade sanitária',v:{p:30,d:5,conf:0.95,N:2500,deff:1,resp:95}},
    {n:'Amostra aleatória simples',v:{p:50,d:5,conf:0.95,N:0,deff:1,resp:100}}
  ],
  refs:['cochran','kish','lohr'],
  compute(v){
    const z=normInv(1-(1-v.conf)/2), p=v.p/100, d=v.d/100;
    if(!(d>0)) return {error:'A margem de erro tem de ser maior que zero.'};
    const n0=z*z*p*(1-p)/(d*d), nF=v.N>0? n0/(1+(n0-1)/v.N):n0, nD=nF*v.deff, nR=nD/(v.resp/100), main=ceil(nR);
    const warn=(d>p||d>1-p)?'A margem de erro é maior do que a própria prevalência: o intervalo incluiria valores abaixo de 0% ou acima de 100%. Reduza a margem.':(v.N>0&&main>v.N?'A amostra final é maior do que a população. Considere um censo.':'');
    const cf=fmt(v.conf*100);
    return {main, mainTxt:fmt(main), label:'entrevistas a planear', sub:`Para estimar ${fmt(v.p)}% ± ${fmt(v.d,1)} p.p. com ${cf}% de confiança`, warn,
      steps:[
        {l:'Amostra base, n₀',n:`z = ${fmt(z,3)} para ${cf}% de confiança`,v:n0},
        {l:'Correcção para população finita',n:v.N>0?`N = ${fmt(v.N)}`:'Não aplicada: população grande',v:nF},
        {l:'× Efeito de desenho',n:`deff = ${fmt(v.deff,1)}`,v:nD},
        {l:'÷ Taxa de resposta',n:`${fmt(v.resp)}% das pessoas respondem`,v:nR,final:true}
      ],
      formula:['n₀ = z² · p(1 − p) / d²','n = n₀ / [1 + (n₀ − 1) / N]','n_final = n · deff / taxa de resposta'],
      fnote:'Cochran (1977), cap. 4; deff segundo Kish (1965).',
      pt:`O tamanho da amostra foi calculado para estimar uma prevalência esperada de ${fmt(v.p)}%, com uma margem de erro de ±${fmt(v.d,1)} pontos percentuais e ${cf}% de confiança, o que resulta numa amostra base de ${ceil(n0)} unidades (Cochran, 1977).${v.N>0?` Aplicou-se a correcção para população finita (N = ${fmt(v.N)}), obtendo-se ${ceil(nF)} unidades.`:''}${v.deff>1?` Para ter em conta o desenho por conglomerados, aplicou-se um efeito de desenho de ${fmt(v.deff,1)} (Kish, 1965).`:''} Considerando uma taxa de resposta esperada de ${fmt(v.resp)}%, a amostra final a planear é de ${fmt(main)} entrevistas.`,
      en:`Sample size was calculated to estimate an expected prevalence of ${fmtEn(v.p)}% with a margin of error of ±${fmtEn(v.d,1)} percentage points at ${fmtEn(v.conf*100)}% confidence, giving a base sample of ${ceil(n0)} units (Cochran, 1977).${v.N>0?` A finite population correction was applied (N = ${fmtEn(v.N)}), giving ${ceil(nF)} units.`:''}${v.deff>1?` To account for the clustered design, a design effect of ${fmtEn(v.deff,1)} was applied (Kish, 1965).`:''} Assuming an expected response rate of ${fmtEn(v.resp)}%, the final planned sample is ${fmtEn(main)} interviews.`};
  },
  sens(v){return {k:'d',xs:range(1,15,0.5),xl:'Margem de erro (p.p.)',yl:'Entrevistas',fx:x=>fmt(x,1),fy:y=>fmt(y)};}
},
{ id:'media', k:'Inquérito', n:'Média', d:'Estimar um valor médio, como rendimento ou pontuação',
  lead:'Para estimar a média de uma variável contínua, como o rendimento mensal do agregado ou a pontuação num teste de leitura.',
  fields:[
    {k:'sigma',l:'Desvio-padrão esperado',u:'unidades',def:15,min:0.01,step:0.5,h:'De um estudo anterior ou piloto. Uma regra prática: amplitude esperada ÷ 4.'},
    {k:'d',l:'Margem de erro (±)',u:'unidades',def:2,min:0.01,step:0.5,h:'Na mesma unidade da variável.'},
    {k:'conf',l:'Nível de confiança',sel:CONF,def:0.95},
    {k:'N',l:'Tamanho da população',u:'pessoas',def:0,min:0,step:100,h:'Deixe 0 se a população for grande.'},
    {k:'deff',l:'Efeito de desenho',u:'deff',def:2,min:1,max:6,step:0.1,range:true},
    {k:'resp',l:'Taxa de resposta esperada',u:'%',def:90,min:30,max:100,step:1,range:true}
  ],
  presets:[
    {n:'Pontuação de leitura (0–100)',v:{sigma:15,d:2,conf:0.95,N:0,deff:2,resp:90}},
    {n:'Rendimento mensal (MT)',v:{sigma:8000,d:1000,conf:0.95,N:0,deff:1.8,resp:85}},
    {n:'Altura de crianças < 5 anos (cm)',v:{sigma:9,d:1,conf:0.95,N:0,deff:1.5,resp:92}}
  ],
  refs:['cochran','kish','lohr'],
  compute(v){
    const z=normInv(1-(1-v.conf)/2);
    if(!(v.d>0&&v.sigma>0)) return {error:'O desvio-padrão e a margem têm de ser maiores que zero.'};
    const n0=z*z*v.sigma*v.sigma/(v.d*v.d), nF=v.N>0?n0/(1+(n0-1)/v.N):n0, nD=nF*v.deff, nR=nD/(v.resp/100), main=ceil(nR);
    const cf=fmt(v.conf*100);
    return {main, mainTxt:fmt(main), label:'entrevistas a planear', sub:`Para estimar a média com ± ${fmt(v.d,2)} unidades e ${cf}% de confiança`, warn:(v.N>0&&main>v.N)?'A amostra final é maior do que a população. Considere um censo.':'',
      steps:[
        {l:'Amostra base, n₀',n:`z = ${fmt(z,3)}; σ = ${fmt(v.sigma,2)}`,v:n0},
        {l:'Correcção para população finita',n:v.N>0?`N = ${fmt(v.N)}`:'Não aplicada: população grande',v:nF},
        {l:'× Efeito de desenho',n:`deff = ${fmt(v.deff,1)}`,v:nD},
        {l:'÷ Taxa de resposta',n:`${fmt(v.resp)}%`,v:nR,final:true}
      ],
      formula:['n₀ = z² · σ² / d²','n = n₀ / [1 + (n₀ − 1) / N]','n_final = n · deff / taxa de resposta'],
      fnote:'Cochran (1977), cap. 2 e 4.',
      pt:`Para estimar a média com uma margem de erro de ±${fmt(v.d,2)} unidades e ${cf}% de confiança, assumindo um desvio-padrão de ${fmt(v.sigma,2)}, a amostra base é de ${ceil(n0)} unidades (Cochran, 1977).${v.N>0?` Com correcção para população finita (N = ${fmt(v.N)}), obtêm-se ${ceil(nF)}.`:''}${v.deff>1?` Aplicou-se um efeito de desenho de ${fmt(v.deff,1)} (Kish, 1965).`:''} Com uma taxa de resposta esperada de ${fmt(v.resp)}%, a amostra final é de ${fmt(main)} entrevistas.`,
      en:`To estimate the mean with a margin of error of ±${fmtEn(v.d,2)} units at ${fmtEn(v.conf*100)}% confidence, assuming a standard deviation of ${fmtEn(v.sigma,2)}, the base sample is ${ceil(n0)} units (Cochran, 1977).${v.N>0?` With a finite population correction (N = ${fmtEn(v.N)}), this becomes ${ceil(nF)}.`:''}${v.deff>1?` A design effect of ${fmtEn(v.deff,1)} was applied (Kish, 1965).`:''} With an expected response rate of ${fmtEn(v.resp)}%, the final sample is ${fmtEn(main)} interviews.`};
  },
  sens(v){const s=v.d/10;return {k:'d',xs:range(+(v.d*0.4).toFixed(4),v.d*2.5,s),xl:'Margem de erro (unidades)',yl:'Entrevistas',fx:x=>fmt(x,2),fy:y=>fmt(y)};}
},
{ id:'duas-proporcoes', k:'Comparação', n:'Duas proporções', d:'Detectar uma diferença entre grupos (sim/não)',
  lead:'Para ensaios ou comparações em que o resultado é binário, como concluir ou não o ano lectivo.',
  fields:[
    {k:'p1',l:'Proporção no grupo de controlo',u:'%',def:40,min:1,max:99,step:1,range:true},
    {k:'p2',l:'Proporção esperada na intervenção',u:'%',def:50,min:1,max:99,step:1,range:true,h:'A diferença entre os dois grupos é o efeito que quer conseguir detectar.'},
    {k:'alpha',l:'Nível de significância (α)',sel:ALPHA,def:0.05},
    {k:'sides',l:'Teste',sel:SIDES,def:2},
    {k:'power',l:'Poder estatístico',sel:POWER,def:0.8},
    {k:'k',l:'Rácio intervenção : controlo',u:'k',def:1,min:0.25,max:4,step:0.25,h:'1 = grupos iguais. 2 = duas pessoas na intervenção por cada uma no controlo.'},
    {k:'deff',l:'Efeito de desenho',u:'deff',def:1,min:1,max:6,step:0.1,range:true,h:'Use acima de 1 se a amostra for recolhida por conglomerados.'},
    {k:'att',l:'Atrito esperado',u:'%',def:10,min:0,max:60,step:1,range:true,h:'Participantes que se perdem até ao inquérito final.'}
  ],
  presets:[
    {n:'Conclusão escolar: 40% → 50%',v:{p1:40,p2:50,alpha:0.05,sides:2,power:0.8,k:1,deff:1,att:10}},
    {n:'Uso de contraceptivos: 25% → 32%',v:{p1:25,p2:32,alpha:0.05,sides:2,power:0.8,k:1,deff:1.5,att:15}},
    {n:'Partos institucionais: 60% → 70%',v:{p1:60,p2:70,alpha:0.05,sides:2,power:0.9,k:1,deff:1.3,att:5}}
  ],
  refs:['fleiss','duflo'],
  compute(v){
    const p1=v.p1/100,p2=v.p2/100,k=v.k;
    if(p1===p2) return {error:'As duas proporções são iguais: não há diferença a detectar.'};
    if(!(k>0)) return {error:'O rácio entre grupos tem de ser maior que zero.'};
    const za=zA(v.alpha,v.sides), zb=zB(v.power), pb=(p1+k*p2)/(1+k), qb=1-pb;
    const n1=Math.pow(za*Math.sqrt((1+1/k)*pb*qb)+zb*Math.sqrt(p1*(1-p1)+p2*(1-p2)/k),2)/Math.pow(p1-p2,2);
    const nD=n1*v.deff, nA=nD/(1-v.att/100), c=ceil(nA), t=ceil(k*nA), main=c+t;
    const diff=Math.abs(v.p2-v.p1), sd=v.sides===2?'bilateral':'unilateral', sdEn=v.sides===2?'two-sided':'one-sided';
    return {main, mainTxt:fmt(main), label:'participantes no total', sub:`Controlo ${fmt(c)} · Intervenção ${fmt(t)} · diferença de ${fmt(diff)} p.p.`, warn:diff<3?'Diferenças abaixo de 3 p.p. exigem amostras muito grandes. Confirme se esse efeito é realista e relevante.':'',
      steps:[
        {l:'Por grupo de controlo (Fleiss)',n:`z_α = ${fmt(za,3)}; z_β = ${fmt(zb,3)}`,v:n1},
        {l:'× Efeito de desenho',n:`deff = ${fmt(v.deff,1)}`,v:nD},
        {l:'÷ (1 − atrito)',n:`${fmt(v.att)}% de perdas`,v:nA},
        {l:'Total: controlo + intervenção',n:`rácio k = ${fmt(k,2)}`,v:main,final:true}
      ],
      formula:['p̄ = (p₁ + k·p₂) / (1 + k)','n₁ = [z_α·√((1 + 1/k)·p̄q̄) + z_β·√(p₁q₁ + p₂q₂/k)]² / (p₁ − p₂)²','n₂ = k·n₁'],
      fnote:'Fleiss, Levin e Paik (2003), cap. 4, sem correcção de continuidade.',
      pt:`Para detectar uma diferença de ${fmt(v.p1)}% no grupo de controlo para ${fmt(v.p2)}% no grupo de intervenção (${fmt(diff)} pontos percentuais), com poder de ${pct(v.power)}% e nível de significância de ${fmt(v.alpha,2)} (teste ${sd}), são necessários ${ceil(n1)} participantes no grupo de controlo e ${ceil(k*n1)} no grupo de intervenção (Fleiss et al., 2003).${v.deff>1?` Aplicou-se um efeito de desenho de ${fmt(v.deff,1)}.`:''} Prevendo um atrito de ${fmt(v.att)}%, a amostra a recrutar é de ${fmt(c)} participantes no controlo e ${fmt(t)} na intervenção, num total de ${fmt(main)}.`,
      en:`To detect a difference from ${fmtEn(v.p1)}% in the control group to ${fmtEn(v.p2)}% in the intervention group (${fmtEn(diff)} percentage points), with ${fmtEn(v.power*100)}% power and a ${sdEn} significance level of ${fmtEn(v.alpha,2)}, ${ceil(n1)} participants are needed in the control group and ${ceil(k*n1)} in the intervention group (Fleiss et al., 2003).${v.deff>1?` A design effect of ${fmtEn(v.deff,1)} was applied.`:''} Allowing for ${fmtEn(v.att)}% attrition, the recruitment target is ${fmtEn(c)} control and ${fmtEn(t)} intervention participants, ${fmtEn(main)} in total.`};
  },
  sens(v){const up=v.p2>v.p1;const xs=up?range(Math.min(v.p1+2,99),Math.min(v.p1+30,99),1):range(Math.max(1,v.p1-30),Math.max(1,v.p1-2),1);return {k:'p2',xs,xl:'Proporção na intervenção (%)',yl:'Participantes',fx:x=>fmt(x),fy:y=>fmt(y)};}
},
{ id:'duas-medias', k:'Comparação', n:'Duas médias', d:'Detectar uma diferença num resultado contínuo',
  lead:'Para ensaios em que o resultado é contínuo, como a pontuação num teste. Inclui o ajuste pela linha de base, que reduz a amostra.',
  fields:[
    {k:'delta',l:'Diferença a detectar',u:'unidades',def:5,min:0.01,step:0.5,h:'Diferença entre as médias dos dois grupos, na unidade do resultado.'},
    {k:'sd',l:'Desvio-padrão do resultado',u:'unidades',def:15,min:0.01,step:0.5},
    {k:'alpha',l:'Nível de significância (α)',sel:ALPHA,def:0.05},
    {k:'sides',l:'Teste',sel:SIDES,def:2},
    {k:'power',l:'Poder estatístico',sel:POWER,def:0.8},
    {k:'k',l:'Rácio intervenção : controlo',u:'k',def:1,min:0.25,max:4,step:0.25},
    {k:'rho',l:'Correlação com a linha de base',u:'ρ',def:0.5,min:0,max:0.95,step:0.05,range:true,h:'Se medir o mesmo resultado antes da intervenção. A amostra cai pelo factor 1 − ρ². Use 0 se não houver linha de base.'},
    {k:'att',l:'Atrito esperado',u:'%',def:10,min:0,max:60,step:1,range:true}
  ],
  presets:[
    {n:'Teste de leitura com linha de base',v:{delta:5,sd:15,alpha:0.05,sides:2,power:0.8,k:1,rho:0.6,att:10}},
    {n:'Peso ao nascer (g)',v:{delta:100,sd:500,alpha:0.05,sides:2,power:0.8,k:1,rho:0,att:5}},
    {n:'Rendimento de negócios de mulheres (MT)',v:{delta:1500,sd:6000,alpha:0.05,sides:2,power:0.8,k:1,rho:0.4,att:15}}
  ],
  refs:['chow','cohen','borm'],
  compute(v){
    if(!(v.delta>0&&v.sd>0)) return {error:'A diferença e o desvio-padrão têm de ser maiores que zero.'};
    if(!(v.k>0)) return {error:'O rácio entre grupos tem de ser maior que zero.'};
    const za=zA(v.alpha,v.sides), zb=zB(v.power);
    const n1=Math.pow(za+zb,2)*v.sd*v.sd*(1+1/v.k)/(v.delta*v.delta), nA=n1*(1-v.rho*v.rho), nT=nA/(1-v.att/100), c=ceil(nT), t=ceil(v.k*nT), main=c+t, dC=v.delta/v.sd;
    const sd=v.sides===2?'bilateral':'unilateral', sdEn=v.sides===2?'two-sided':'one-sided';
    return {main, mainTxt:fmt(main), label:'participantes no total', sub:`Controlo ${fmt(c)} · Intervenção ${fmt(t)} · d de Cohen = ${fmt(dC,2)}`, warn:dC<0.1?'Um efeito abaixo de 0,1 desvios-padrão exige amostras muito grandes.':'',
      steps:[
        {l:'Por grupo de controlo, sem ajuste',n:`z_α = ${fmt(za,3)}; z_β = ${fmt(zb,3)}`,v:n1},
        {l:'× (1 − ρ²): ajuste pela linha de base',n:v.rho>0?`ρ = ${fmt(v.rho,2)}: poupa ${fmt(v.rho*v.rho*100)}%`:'Sem linha de base',v:nA},
        {l:'÷ (1 − atrito)',n:`${fmt(v.att)}% de perdas`,v:nT},
        {l:'Total: controlo + intervenção',n:`rácio k = ${fmt(v.k,2)}`,v:main,final:true}
      ],
      formula:['n₁ = (z_α + z_β)² · σ² · (1 + 1/k) / δ²','n₁,ANCOVA = n₁ · (1 − ρ²)','n₂ = k·n₁'],
      fnote:'Chow, Shao e Wang (2008), cap. 3; ajuste de Borm et al. (2007).',
      pt:`Para detectar uma diferença de ${fmt(v.delta,2)} unidades entre grupos (d de Cohen = ${fmt(dC,2)}), assumindo um desvio-padrão de ${fmt(v.sd,2)}, com poder de ${pct(v.power)}% e nível de significância de ${fmt(v.alpha,2)} (teste ${sd}), são necessários ${ceil(n1)} participantes no grupo de controlo (Chow et al., 2008).${v.rho>0?` A análise ajustará pela medição de linha de base (ANCOVA); com uma correlação esperada de ${fmt(v.rho,2)}, a amostra reduz-se pelo factor 1 − ρ² para ${ceil(nA)} por grupo (Borm et al., 2007).`:''} Prevendo um atrito de ${fmt(v.att)}%, a amostra a recrutar é de ${fmt(c)} participantes no controlo e ${fmt(t)} na intervenção, num total de ${fmt(main)}.`,
      en:`To detect a between-group difference of ${fmtEn(v.delta,2)} units (Cohen's d = ${fmtEn(dC,2)}), assuming a standard deviation of ${fmtEn(v.sd,2)}, with ${fmtEn(v.power*100)}% power and a ${sdEn} significance level of ${fmtEn(v.alpha,2)}, ${ceil(n1)} participants are needed in the control group (Chow et al., 2008).${v.rho>0?` The analysis will adjust for the baseline measurement (ANCOVA); with an expected correlation of ${fmtEn(v.rho,2)}, the sample is reduced by the factor 1 − ρ² to ${ceil(nA)} per group (Borm et al., 2007).`:''} Allowing for ${fmtEn(v.att)}% attrition, the recruitment target is ${fmtEn(c)} control and ${fmtEn(t)} intervention participants, ${fmtEn(main)} in total.`};
  },
  sens(v){return {k:'delta',xs:range(+(v.delta*0.4).toFixed(4),v.delta*2.5,v.delta/10),xl:'Diferença a detectar (unidades)',yl:'Participantes',fx:x=>fmt(x,1),fy:y=>fmt(y)};}
},
{ id:'conglomerados', k:'Avaliação de impacto', n:'RCT por conglomerados', d:'Escolas, aldeias ou unidades sanitárias aleatorizadas',
  lead:'Quando se aleatorizam grupos inteiros, como escolas, a semelhança entre alunos da mesma escola (ICC) aumenta a amostra necessária.',
  fields:[
    {k:'d',l:'Efeito mínimo a detectar (MDES)',u:'DP',def:0.2,min:0.02,max:1.5,step:0.01,range:true,h:'Em desvios-padrão. Em educação, 0,1 a 0,2 DP já é um efeito relevante.'},
    {k:'icc',l:'Correlação intraclasse (ICC)',u:'ρ',def:0.2,min:0,max:0.6,step:0.01,range:true,h:'Aprendizagem em escolas: 0,15 a 0,35. Saúde em comunidades: 0,01 a 0,10.'},
    {k:'m',l:'Participantes por conglomerado',u:'por escola',def:20,min:2,step:1},
    {k:'r2',l:'Variância explicada por covariáveis',u:'R²',def:0.4,min:0,max:0.9,step:0.05,range:true,h:'Com um teste de linha de base, R² de 0,3 a 0,6 é comum. Use 0 se não houver.'},
    {k:'alpha',l:'Nível de significância (α, bilateral)',sel:ALPHA,def:0.05},
    {k:'power',l:'Poder estatístico',sel:POWER,def:0.8},
    {k:'att',l:'Atrito dentro do conglomerado',u:'%',def:10,min:0,max:60,step:1,range:true}
  ],
  presets:[
    {n:'Aprendizagem em escolas primárias',v:{d:0.2,icc:0.2,m:20,r2:0.4,alpha:0.05,power:0.8,att:10}},
    {n:'Saúde comunitária (aldeias)',v:{d:0.25,icc:0.05,m:25,r2:0,alpha:0.05,power:0.8,att:10}},
    {n:'Grupos de poupança de mulheres',v:{d:0.25,icc:0.1,m:15,r2:0.2,alpha:0.05,power:0.8,att:15}}
  ],
  refs:['hayes','donner','bloom05','dong','duflo'],
  compute(v){
    const mE=v.m*(1-v.att/100);
    if(!(v.d>0)) return {error:'O efeito a detectar tem de ser maior que zero.'};
    if(mE<1) return {error:'Com este atrito sobra menos de um participante por conglomerado.'};
    const P=0.5, za=zA(v.alpha,2), zb=zB(v.power);
    const base=(v.icc+(1-v.icc)/mE)*(1-v.r2)/(P*(1-P)*v.d*v.d);
    const nInd=Math.pow(za+zb,2)*(1-v.r2)/(P*(1-P)*v.d*v.d), deff=1+(mE-1)*v.icc;
    const J0=Math.pow(za+zb,2)*base; let J=J0, df=1;
    for(let i=0;i<8;i++){df=Math.max(J-2,2);const M=tInv(1-v.alpha/2,df)+tInv(v.power,df);J=M*M*base;}
    const arm=Math.max(2,ceil(J/2)), tot=arm*2, people=tot*v.m;
    return {main:arm, mainTxt:fmt(arm), label:'conglomerados por braço', sub:`${fmt(tot)} conglomerados no total · ${fmt(people)} participantes a recrutar`, warn:tot<10?'Com menos de 10 conglomerados, a aleatorização pode não equilibrar os grupos. Considere estratificar ou usar mais conglomerados.':(v.icc>=0.3?'Com um ICC tão alto, acrescentar conglomerados ajuda muito mais do que acrescentar participantes por conglomerado.':''),
      steps:[
        {l:'Amostra total se a aleatorização fosse individual',n:`já com R² = ${fmt(v.r2,2)}`,v:nInd},
        {l:'× Efeito de desenho',n:`1 + (${fmt(mE,1)} − 1) × ${fmt(v.icc,2)} = ${fmt(deff,2)}`,v:nInd*deff},
        {l:'÷ participantes efectivos por conglomerado',n:`${fmt(v.m)} × (1 − ${fmt(v.att)}%) = ${fmt(mE,1)}; conglomerados (aprox. normal)`,v:J0,bar:false},
        {l:'Total de conglomerados, com correcção t',n:`gl = J − 2 ≈ ${fmt(df,0)}; arredondado para ${fmt(tot)}`,v:tot,bar:false,final:true}
      ],
      formula:['deff = 1 + (m − 1)·ICC','J = M² · [ICC + (1 − ICC)/m] · (1 − R²) / [P(1 − P)·MDES²]','M = t(1 − α/2; J − 2) + t(1 − β; J − 2)'],
      fnote:'Bloom (2005); Dong e Maynard (2013); Hayes e Moulton (2017). P = 0,5 (braços iguais).',
      pt:`O estudo é um ensaio aleatorizado por conglomerados com dois braços de igual dimensão. Para detectar um efeito mínimo de ${fmt(v.d,2)} desvios-padrão, com poder de ${pct(v.power)}% e nível de significância de ${fmt(v.alpha,2)} (bilateral), assumiu-se um coeficiente de correlação intraclasse de ${fmt(v.icc,2)}, ${fmt(v.m)} participantes por conglomerado, ${fmt(v.att)}% de atrito e covariáveis de linha de base que explicam ${fmt(v.r2*100)}% da variância. O efeito de desenho correspondente é de ${fmt(deff,2)} (Donner e Klar, 2000). São necessários ${fmt(arm)} conglomerados por braço (${fmt(tot)} no total; ${fmt(people)} participantes). O cálculo segue Bloom (2005) e Dong e Maynard (2013), com multiplicador t para ${fmt(Math.max(tot-2,1))} graus de liberdade.`,
      en:`The study is a two-arm cluster randomised trial with equal allocation. To detect a minimum effect of ${fmtEn(v.d,2)} standard deviations with ${fmtEn(v.power*100)}% power and a two-sided significance level of ${fmtEn(v.alpha,2)}, we assumed an intraclass correlation of ${fmtEn(v.icc,2)}, ${fmtEn(v.m)} participants per cluster, ${fmtEn(v.att)}% attrition and baseline covariates explaining ${fmtEn(v.r2*100)}% of the variance. The corresponding design effect is ${fmtEn(deff,2)} (Donner and Klar, 2000). ${fmtEn(arm)} clusters per arm are required (${fmtEn(tot)} in total; ${fmtEn(people)} participants). Calculations follow Bloom (2005) and Dong and Maynard (2013), using a t multiplier with ${fmtEn(Math.max(tot-2,1))} degrees of freedom.`};
  },
  sens(v){return {k:'icc',xs:range(0.01,0.4,0.01),xl:'ICC',yl:'Conglomerados por braço',fx:x=>fmt(x,2),fy:y=>fmt(y)};}
},
{ id:'mdes', k:'Modo inverso', n:'Efeito mínimo detectável', d:'Tenho N conglomerados. O que consigo detectar?',
  lead:'Quando o orçamento já fixa o número de escolas ou aldeias, este módulo diz qual o menor efeito que o estudo consegue detectar.',
  fields:[
    {k:'J',l:'Total de conglomerados',u:'escolas',def:60,min:4,step:2},
    {k:'m',l:'Participantes por conglomerado',u:'por escola',def:20,min:2,step:1},
    {k:'icc',l:'Correlação intraclasse (ICC)',u:'ρ',def:0.2,min:0,max:0.6,step:0.01,range:true},
    {k:'r2',l:'Variância explicada por covariáveis',u:'R²',def:0.4,min:0,max:0.9,step:0.05,range:true},
    {k:'P',l:'Conglomerados no tratamento',u:'%',def:50,min:10,max:90,step:5,range:true},
    {k:'alpha',l:'Nível de significância (α, bilateral)',sel:ALPHA,def:0.05},
    {k:'power',l:'Poder estatístico',sel:POWER,def:0.8},
    {k:'att',l:'Atrito dentro do conglomerado',u:'%',def:10,min:0,max:60,step:1,range:true}
  ],
  presets:[
    {n:'60 escolas, 20 alunas por escola',v:{J:60,m:20,icc:0.2,r2:0.4,P:50,alpha:0.05,power:0.8,att:10}},
    {n:'30 aldeias, 25 agregados',v:{J:30,m:25,icc:0.05,r2:0,P:50,alpha:0.05,power:0.8,att:10}},
    {n:'100 escolas, 3 braços (1/3 tratamento)',v:{J:100,m:15,icc:0.2,r2:0.4,P:35,alpha:0.05,power:0.8,att:10}}
  ],
  refs:['bloom95','bloom05','dong','kraft'],
  compute(v){
    const mE=v.m*(1-v.att/100), P=v.P/100;
    if(!(v.J>=4)) return {error:'São precisos pelo menos 4 conglomerados.'};
    if(mE<1) return {error:'Com este atrito sobra menos de um participante por conglomerado.'};
    const df=v.J-2, M=tInv(1-v.alpha/2,df)+tInv(v.power,df), comp=(v.icc+(1-v.icc)/mE)*(1-v.r2), se=Math.sqrt(comp/(P*(1-P)*v.J)), mdes=M*se;
    const interp=mdes<0.05?'Muito sensível: detecta efeitos pequenos segundo Kraft (2020).':mdes<0.2?'Detecta efeitos médios em educação (0,05 a 0,20 DP segundo Kraft, 2020). Boa sensibilidade.':mdes<0.35?'Só detecta efeitos grandes (≥ 0,20 DP). Muitas intervenções eficazes em educação ficariam abaixo deste limiar.':'Pouca sensibilidade: o estudo só detectaria efeitos muito grandes. Aumente o número de conglomerados.';
    return {main:mdes, mainTxt:fmt(mdes,2), label:'desvios-padrão: efeito mínimo detectável', sub:interp, warn:v.J<10?'Com menos de 10 conglomerados as estimativas de poder são pouco fiáveis.':'',
      steps:[
        {l:'Multiplicador M',n:`t(${fmt(1-v.alpha/2,3)}; ${df}) + t(${fmt(v.power,2)}; ${df}) graus de liberdade`,v:M,dec:3,bar:false},
        {l:'Variância por conglomerado',n:`[ICC + (1 − ICC)/m] × (1 − R²), com m = ${fmt(mE,1)}`,v:comp,dec:4,bar:false},
        {l:'Erro-padrão do efeito',n:`√(variância / [P(1 − P)·J])`,v:se,dec:4,bar:false},
        {l:'MDES = M × erro-padrão',n:'em desvios-padrão do resultado',v:mdes,dec:3,bar:false,final:true}
      ],
      formula:['MDES = M · √{ [ICC + (1 − ICC)/m] · (1 − R²) / [P(1 − P)·J] }','M = t(1 − α/2; J − 2) + t(1 − β; J − 2)'],
      fnote:'Bloom (1995, 2005); Dong e Maynard (2013). Com J grande, M ≈ 2,8 para 80% de poder.',
      pt:`Com ${fmt(v.J)} conglomerados (${fmt(v.P)}% atribuídos ao tratamento) e ${fmt(v.m)} participantes por conglomerado, assumindo um coeficiente de correlação intraclasse de ${fmt(v.icc,2)}, ${fmt(v.att)}% de atrito e covariáveis que explicam ${fmt(v.r2*100)}% da variância, o efeito mínimo detectável é de ${fmt(mdes,2)} desvios-padrão, com poder de ${pct(v.power)}% e nível de significância de ${fmt(v.alpha,2)} (bilateral) (Bloom, 1995; Dong e Maynard, 2013).`,
      en:`With ${fmtEn(v.J)} clusters (${fmtEn(v.P)}% assigned to treatment) and ${fmtEn(v.m)} participants per cluster, assuming an intraclass correlation of ${fmtEn(v.icc,2)}, ${fmtEn(v.att)}% attrition and covariates explaining ${fmtEn(v.r2*100)}% of the variance, the minimum detectable effect size is ${fmtEn(mdes,2)} standard deviations at ${fmtEn(v.power*100)}% power and a two-sided significance level of ${fmtEn(v.alpha,2)} (Bloom, 1995; Dong and Maynard, 2013).`};
  },
  sens(v){return {k:'J',xs:range(10,200,10),xl:'Total de conglomerados',yl:'MDES (DP)',fx:x=>fmt(x),fy:y=>fmt(y,2)};}
}
];

/* ---------- state ---------- */
const state={}; MODS.forEach(m=>{state[m.id]={};m.fields.forEach(f=>state[m.id][f.k]=f.def);});
let cur=MODS[0].id, lang='pt';
try{const h=location.hash.slice(1); if(MODS.some(m=>m.id===h)) cur=h;}catch(e){}
const $=id=>document.getElementById(id);
const esc=s=>String(s).replace(/[&<>"]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));

function renderRail(){
  $('rail').innerHTML=MODS.map(m=>`<button class="mod" type="button" data-id="${m.id}" aria-pressed="${m.id===cur}"><span class="k">${m.k}</span><span class="n">${m.n}</span><span class="d">${m.d}</span></button>`).join('');
  $('rail').querySelectorAll('.mod').forEach(b=>b.addEventListener('click',()=>{cur=b.dataset.id;try{history.replaceState(null,'','#'+cur);}catch(e){}renderRail();renderInputs();recalc();}));
}
function fieldHTML(m,f){
  const id=`f-${m.id}-${f.k}`, val=state[m.id][f.k];
  let ctl;
  if(f.sel){ctl=`<select id="${id}" data-k="${f.k}">${f.sel.map(([v,t])=>`<option value="${v}"${+v===+val?' selected':''}>${t}</option>`).join('')}</select>`;}
  else if(f.range){ctl=`<input type="range" id="${id}-r" data-k="${f.k}" min="${f.min}" max="${f.max}" step="${f.step}" value="${val}" aria-label="${esc(f.l)}"><input type="number" class="sm" id="${id}" data-k="${f.k}" min="${f.min}" max="${f.max}" step="${f.step}" value="${val}" inputmode="decimal">`;}
  else{ctl=`<input type="number" id="${id}" data-k="${f.k}" ${f.min!=null?`min="${f.min}"`:''} ${f.max!=null?`max="${f.max}"`:''} step="${f.step||'any'}" value="${val}" inputmode="decimal">`;}
  return `<div class="field"><label for="${id}"><span>${f.l}</span>${f.u?`<span class="u">${f.u}</span>`:''}</label><div class="ctl">${ctl}</div>${f.h?`<div class="h">${f.h}</div>`:''}</div>`;
}
function renderInputs(){
  const m=MODS.find(x=>x.id===cur);
  $('inputs').innerHTML=`<h2>${m.n}</h2><p class="lead">${m.lead}</p>
   <div class="presets"><span>Exemplos (valores indicativos)</span>${m.presets.map((p,i)=>`<button type="button" class="chip" data-p="${i}">${p.n}</button>`).join('')}</div>
   <div class="fields">${m.fields.map(f=>fieldHTML(m,f)).join('')}</div>`;
  $('inputs').querySelectorAll('[data-k]').forEach(el=>el.addEventListener('input',()=>{
    const k=el.dataset.k, x=parseFloat(el.value); state[cur][k]=x;
    const f=m.fields.find(q=>q.k===k);
    if(f.range){const twin=el.type==='range'?$(`f-${m.id}-${k}`):$(`f-${m.id}-${k}-r`); if(twin&&!isNaN(x)) twin.value=x;}
    recalc();
  }));
  $('inputs').querySelectorAll('[data-p]').forEach(b=>b.addEventListener('click',()=>{Object.assign(state[cur],m.presets[+b.dataset.p].v);renderInputs();recalc();}));
}

/* ---------- chart ---------- */
function niceMax(v){if(!(v>0))return 1;const e=Math.pow(10,Math.floor(Math.log10(v))),f=v/e;return (f<=1?1:f<=2?2:f<=2.5?2.5:f<=5?5:10)*e;}
function chart(pts,cx,cy,s){
  const W=560,H=240,L=58,R=20,T=20,B=42;
  const xs=pts.map(p=>p[0]), xmin=Math.min(...xs), xmax=Math.max(...xs), ymax=niceMax(Math.max(...pts.map(p=>p[1]),cy));
  const X=x=>L+(x-xmin)/((xmax-xmin)||1)*(W-L-R), Y=y=>T+(1-y/ymax)*(H-T-B);
  let g='';
  for(let i=0;i<=4;i++){const yv=ymax*i/4,yy=Y(yv);g+=`<line class="g" x1="${L}" x2="${W-R}" y1="${yy}" y2="${yy}"/><text x="${L-8}" y="${yy+4}" text-anchor="end">${s.fy(yv)}</text>`;}
  for(let i=0;i<=4;i++){const xv=xmin+(xmax-xmin)*i/4,xx=X(xv);g+=`<text x="${xx}" y="${H-B+18}" text-anchor="middle">${s.fx(xv)}</text>`;}
  const d=pts.map((p,i)=>`${i?'L':'M'}${X(p[0]).toFixed(1)},${Y(Math.min(p[1],ymax)).toFixed(1)}`).join('');
  const area=`${d}L${X(xmax).toFixed(1)},${Y(0)}L${X(xmin).toFixed(1)},${Y(0)}Z`;
  const px=X(cx),py=Y(cy),anchor=px>W-130?'end':'start',dx=anchor==='end'?-10:10;
  return `<svg class="chart" viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(s.yl)} em função de ${esc(s.xl)}">
   ${g}<path class="ar" d="${area}"/><path class="ln" d="${d}"/>
   <line class="g" x1="${px}" x2="${px}" y1="${py}" y2="${Y(0)}" stroke-dasharray="3 3"/>
   <circle class="dot" cx="${px}" cy="${py}" r="5.5"/>
   <text class="dl" x="${px+dx}" y="${Math.max(py-10,T+10)}" text-anchor="${anchor}">Actual: ${s.fy(cy)}</text>
   <text x="${(L+W-R)/2}" y="${H-6}" text-anchor="middle">${esc(s.xl)}</text>
  </svg>`;
}

/* ---------- recalc ---------- */
function recalc(){
  const m=MODS.find(x=>x.id===cur), v=state[cur];
  if(m.fields.some(f=>isNaN(v[f.k]))){$('results').innerHTML=`<div class="warn">Preencha todos os campos com números para ver o resultado.</div>`;return;}
  const r=m.compute(v);
  if(r.error){$('results').innerHTML=`<div class="warn">${r.error}</div>`;return;}
  const s=m.sens(v); let xs=s.xs.slice(); const cx=v[s.k]; if(!xs.some(x=>Math.abs(x-cx)<1e-9)) xs.push(cx); xs.sort((a,b)=>a-b);
  const pts=xs.map(x=>{const q=m.compute({...v,[s.k]:x});return q.error?null:[x,q.main];}).filter(Boolean);
  const barMax=Math.max(...r.steps.filter(st=>st.bar!==false).map(st=>st.v),1);
  const steps=r.steps.map(st=>`<div class="step${st.final?' final':''}"><span class="sl">${st.l}</span><span class="sv">${st.dec?fmt(st.v,st.dec):fmt(st.v,st.v<100&&!st.final?1:0)}</span><span class="sn">${st.n}</span>${st.bar===false?'':`<span class="bar"><i style="width:${Math.min(100,st.v/barMax*100).toFixed(1)}%"></i></span>`}</div>`).join('');
  $('results').innerHTML=`
    <div class="big"><span class="v">${r.mainTxt}</span><div><div class="l">${r.label}</div><div class="s">${r.sub}</div></div></div>
    ${r.warn?`<div class="warn">${r.warn}</div>`:''}
    <div><div class="sec-t">Do cálculo ao terreno</div><div class="steps">${steps}</div></div>
    <div><div class="sec-t">Fórmula</div><div class="formula">${r.formula.map(f=>`<div>${esc(f)}</div>`).join('')}<div class="c">${r.fnote}</div></div></div>
    <div><div class="sec-t">Sensibilidade<span style="text-transform:none;letter-spacing:0;font-weight:400">${esc(s.yl)} conforme ${esc(s.xl.toLowerCase())}</span></div>${pts.length>1?chart(pts,cx,r.main,s):''}</div>
    <div><div class="sec-t">Texto para o protocolo
      <span style="display:flex;gap:8px;align-items:center"><span class="seg" role="group" aria-label="Idioma"><button type="button" data-l="pt" aria-pressed="${lang==='pt'}">PT</button><button type="button" data-l="en" aria-pressed="${lang==='en'}">EN</button></span><button type="button" class="btn" id="copy">Copiar</button></span></div>
      <p class="proto" id="proto">${esc(lang==='pt'?r.pt:r.en)}</p></div>
    <div><div class="sec-t">Obras que inspiraram este cálculo</div><ul class="refs">${m.refs.map(k=>`<li>${REFS[k].f}</li>`).join('')}</ul></div>`;
  $('results').querySelectorAll('[data-l]').forEach(b=>b.addEventListener('click',()=>{lang=b.dataset.l;recalc();}));
  $('copy').addEventListener('click',()=>{
    const t=$('proto').textContent, btn=$('copy');
    const sel=()=>{const rg=document.createRange();rg.selectNodeContents($('proto'));const ss=getSelection();ss.removeAllRanges();ss.addRange(rg);btn.textContent='Seleccionado: use Ctrl+C';};
    try{navigator.clipboard.writeText(t).then(()=>{btn.textContent='Copiado';setTimeout(()=>{btn.textContent='Copiar';},1800);}).catch(sel);}catch(e){sel();}
  });
}

/* ---------- library ---------- */
$('libtable').innerHTML=Object.values(REFS).map(r=>`<tr><td>${r.f}</td><td>${r.u}</td><td>${r.m}</td></tr>`).join('');

renderRail(); renderInputs(); recalc();

/* ================= v2: references for new tools ================= */
Object.assign(REFS,{
  neyman:{c:'Neyman, 1934',f:'Neyman, J. (1934). On the two different aspects of the representative method. <i>Journal of the Royal Statistical Society</i>, 97(4), 558–625.',u:'Alocação óptima entre estratos.',m:'Estratificação'},
  glenn:{c:'Glennerster e Takavarasha, 2013',f:'Glennerster, R. e Takavarasha, K. (2013). <i>Running Randomized Evaluations: A Practical Guide</i>. Princeton University Press.',u:'Aleatorização estratificada e prática de terreno.',m:'Aleatorização'},
  bruhn:{c:'Bruhn e McKenzie, 2009',f:'Bruhn, M. e McKenzie, D. (2009). In pursuit of balance: Randomization in practice in development field experiments. <i>American Economic Journal: Applied Economics</i>, 1(4), 200–232.',u:'Equilíbrio entre grupos através de sorteio dentro de estratos.',m:'Aleatorização'}
});
REFS.cochran.m+=', Estratificação'; REFS.lohr.m+=', Estratificação'; REFS.kish.m+=', Estratificação';
$('libtable').innerHTML=Object.values(REFS).map(r=>`<tr><td>${r.f}</td><td>${r.u}</td><td>${r.m}</td></tr>`).join('');

const norm=s=>String(s).normalize('NFD').replace(/[̀-ͯ]/g,'').toLowerCase();
function download(name,text){
  try{const blob=new Blob(['﻿'+text],{type:'text/csv;charset=utf-8'});const a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download=name;document.body.appendChild(a);a.click();setTimeout(()=>{URL.revokeObjectURL(a.href);a.remove();},500);}catch(e){}
}
const csvCell=v=>{const s=String(v);return /[;"\n]/.test(s)?'"'+s.replace(/"/g,'""')+'"':s;};
function copyText(btn,getText,fallbackEl){
  const sel=()=>{const rg=document.createRange();rg.selectNodeContents(fallbackEl);const ss=getSelection();ss.removeAllRanges();ss.addRange(rg);btn.textContent='Seleccionado: use Ctrl+C';};
  try{navigator.clipboard.writeText(getText()).then(()=>{btn.textContent='Copiado';setTimeout(()=>{btn.textContent='Copiar';},1800);}).catch(sel);}catch(e){sel();}
}

/* ================= data ================= */
let DATA=null;
fetch('data/mz-censo-2017.json').then(r=>r.json()).then(d=>{
  DATA=d; DATA.national=d.provinces.reduce((s,p)=>s+p.pop,0);
  DATA.provPop={}; d.provinces.forEach(p=>DATA.provPop[p.n]=p.pop);
  $('src-line').innerHTML=`${esc(d.source)}. Data do censo: 1 de Agosto de 2017. ${esc(d.note)} <a href="${d.sourceUrl}" target="_blank" rel="noopener">Ver fonte</a>.`;
  $('q-source').innerHTML=`Fonte: ${esc(d.source)}. ${esc(d.note)}`;
  initStrata(); initSearch();
}).catch(()=>{$('st-results').innerHTML='<div class="warn">Não foi possível carregar os dados do Censo 2017. Recarregue a página.</div>';});

/* ================= stratification ================= */
const ST={level:'prov',prov:'Nampula',sel:new Set(),pH:{},p:50,d:5,conf:0.95,deff:2,resp:90,min:30,target:'global',alloc:'prop',cost:1200,perday:6,team:12,filter:'',lang:'pt'};
const keyOf=u=>u.kind+':'+u.p+'|'+u.n;
function unitsInView(){
  if(ST.level==='prov') return DATA.provinces.map(p=>({kind:'P',p:p.n,n:p.n,pop:p.pop}));
  return DATA.districts.filter(x=>ST.prov==='*'||x.p===ST.prov).map(x=>({kind:'D',p:x.p,n:x.n,pop:x.pop,city:x.city,from:x.from}));
}
function selectAllInView(){ST.sel=new Set(unitsInView().map(keyOf));}
const ALLOC_H={prop:'Cada estrato recebe entrevistas na proporção da população. Óptima para o total nacional.',raiz:'Proporcional à raiz quadrada da população (Kish). Reforça os estratos pequenos com pouca perda no total.',igual:'O mesmo número em todos os estratos. Bom para comparar estratos entre si.',neyman:'Proporcional à população × desvio-padrão do estrato (Neyman). Edite a prevalência esperada de cada estrato na tabela.'};
const TARGET_H={global:'Calcula a amostra para a margem pedida no total e distribui-a pelos estratos.',estrato:'Cada estrato recebe a amostra necessária para ter, sozinho, a margem pedida. A amostra total cresce muito.'};

function initStrata(){
  $('st-prov').innerHTML=`<option value="*">Todas as províncias (161 unidades)</option>`+DATA.provinces.map(p=>`<option value="${esc(p.n)}">${esc(p.n)}</option>`).join('');
  [['st-p','p'],['st-d','d'],['st-deff','deff'],['st-resp','resp'],['st-min','min'],['st-cost','cost'],['st-perday','perday'],['st-team','team']].forEach(([id,k])=>{const el=$(id);el.value=ST[k];el.addEventListener('input',()=>{ST[k]=parseFloat(el.value);calcStrata();});});
  $('st-conf').value=String(ST.conf); $('st-conf').addEventListener('change',e=>{ST.conf=+e.target.value;calcStrata();});
  $('st-prov').value=ST.prov; $('st-prov').addEventListener('change',e=>{ST.prov=e.target.value;selectAllInView();renderUnits();calcStrata();});
  const seg=(id,k,after)=>{$(id).querySelectorAll('button').forEach(b=>b.addEventListener('click',()=>{ST[k]=b.dataset.v;if(after)after();syncSeg();calcStrata();}));};
  seg('st-level','level',()=>{selectAllInView();renderUnits();});
  seg('st-target','target'); seg('st-alloc','alloc');
  $('st-search').addEventListener('input',e=>{ST.filter=e.target.value;renderUnits();});
  $('st-all').addEventListener('click',()=>{filteredUnits().forEach(u=>ST.sel.add(keyOf(u)));renderUnits();calcStrata();});
  $('st-none').addEventListener('click',()=>{filteredUnits().forEach(u=>ST.sel.delete(keyOf(u)));renderUnits();calcStrata();});
  selectAllInView(); syncSeg(); renderUnits(); calcStrata();
}
function syncSeg(){
  [['st-level','level'],['st-target','target'],['st-alloc','alloc']].forEach(([id,k])=>$(id).querySelectorAll('button').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.v===ST[k]))));
  $('st-prov-wrap').hidden=ST.level!=='dist'; $('st-prov').value=ST.prov;
  $('st-alloc-h').textContent=ALLOC_H[ST.alloc]; $('st-target-h').textContent=TARGET_H[ST.target];
  $('st-alloc').style.opacity=ST.target==='estrato'?.45:1;
}
function filteredUnits(){const f=norm(ST.filter.trim());return unitsInView().filter(u=>!f||norm(u.n).includes(f)||norm(u.p).includes(f));}
function renderUnits(){
  const us=filteredUnits(), all=unitsInView();
  $('st-units').innerHTML=us.length?us.map(u=>{const k=keyOf(u);return `<label><span><input type="checkbox" data-key="${esc(k)}"${ST.sel.has(k)?' checked':''}>${esc(u.n)}${ST.level==='dist'&&ST.prov==='*'?`<span class="tagc">${esc(u.p)}</span>`:''}</span><span class="pp">${fmt(u.pop)}</span></label>`;}).join(''):'<p class="note" style="padding:8px 10px;margin:0">Nenhuma unidade com esse nome.</p>';
  $('st-units').querySelectorAll('input').forEach(c=>c.addEventListener('change',()=>{c.checked?ST.sel.add(c.dataset.key):ST.sel.delete(c.dataset.key);updateCount();calcStrata();}));
  updateCount();
  function updateCount(){const n=all.filter(u=>ST.sel.has(keyOf(u))).length;$('st-count').textContent=`${n} de ${all.length} seleccionados`;}
}

function stratCompute(){
  const units=unitsInView().filter(u=>ST.sel.has(keyOf(u)));
  if(!units.length) return {error:'Seleccione pelo menos um estrato.'};
  if([ST.p,ST.d,ST.deff,ST.resp,ST.min].some(x=>isNaN(x))||ST.d<=0) return {error:'Preencha os parâmetros com números válidos.'};
  const z=normInv(1-(1-ST.conf)/2), d=ST.d/100, N=units.reduce((s,u)=>s+u.pop,0), resp=ST.resp/100;
  const rows=units.map(u=>{const k=keyOf(u);const ph=(ST.pH[k]??ST.p)/100;return {...u,key:k,ph,W:u.pop/N,S:Math.sqrt(ph*(1-ph))};});
  let nTot=0;
  if(ST.target==='global'){
    const n0=z*z*(ST.p/100)*(1-ST.p/100)/(d*d)*ST.deff, nG=n0/(1+(n0-1)/N);
    const w=r=>ST.alloc==='prop'?r.pop:ST.alloc==='raiz'?Math.sqrt(r.pop):ST.alloc==='igual'?1:r.pop*r.S;
    const sw=rows.reduce((s,r)=>s+w(r),0);
    rows.forEach(r=>{r.nRaw=nG*w(r)/sw;}); nTot=nG;
  } else {
    rows.forEach(r=>{const n0=z*z*r.ph*(1-r.ph)/(d*d)*ST.deff;r.nRaw=n0/(1+(n0-1)/r.pop);});
  }
  rows.forEach(r=>{r.ns=Math.min(r.pop,Math.max(ST.min,Math.ceil(r.nRaw)));r.nField=Math.ceil(r.ns/resp);r.minHit=r.ns>Math.ceil(r.nRaw);
    const fpc=Math.max(0,1-r.ns/r.pop); r.marg=z*Math.sqrt(r.ph*(1-r.ph)*ST.deff/r.ns*fpc);});
  const V=rows.reduce((s,r)=>s+r.W*r.W*r.ph*(1-r.ph)*ST.deff/r.ns*Math.max(0,1-r.ns/r.pop),0);
  const overall=z*Math.sqrt(V), net=rows.reduce((s,r)=>s+r.ns,0), field=rows.reduce((s,r)=>s+r.nField,0);
  const cost=field*(ST.cost||0), enumDays=field/(ST.perday||1), days=Math.ceil(enumDays/(ST.team||1));
  return {rows,N,z,net,field,overall,cost,enumDays,days,minHits:rows.filter(r=>r.minHit).length,nTot};
}
const ALLOC_PT={prop:'proporcional à população',raiz:'proporcional à raiz quadrada da população (Kish, 1965)',igual:'igual entre estratos',neyman:'óptima de Neyman (Neyman, 1934)'};
const ALLOC_EN={prop:'proportional to population',raiz:'proportional to the square root of population (Kish, 1965)',igual:'equal across strata',neyman:'Neyman optimal allocation (Neyman, 1934)'};
function calcStrata(){
  if(!DATA) return;
  const R=stratCompute();
  if(R.error){$('st-results').innerHTML=`<div class="warn">${R.error}</div>`;return;}
  const lvl=ST.level==='prov'?'províncias':'distritos', lvlEn=ST.level==='prov'?'provinces':'districts';
  const scope=ST.level==='dist'&&ST.prov!=='*'?` da província de ${ST.prov}`:'', scopeEn=ST.level==='dist'&&ST.prov!=='*'?` of ${ST.prov} province`:'';
  const H=R.rows.length, bad=R.rows.filter(r=>r.marg>0.10).length, maxM=Math.max(...R.rows.map(r=>r.marg));
  const cf=fmt(ST.conf*100);
  const pt=`A amostra foi estratificada por ${H} ${lvl}${scope}, com base na população do Censo 2017 (INE). ${ST.target==='global'?`Para estimar uma prevalência de ${fmt(ST.p)}% com margem de erro de ±${fmt(ST.d,1)} pontos percentuais e ${cf}% de confiança, com um efeito de desenho de ${fmt(ST.deff,1)}, calculou-se uma amostra total de ${fmt(Math.ceil(R.nTot))} entrevistas (Cochran, 1977), distribuída de forma ${ALLOC_PT[ST.alloc]}`:`Para que cada estrato tenha uma margem de erro de ±${fmt(ST.d,1)} pontos percentuais com ${cf}% de confiança e um efeito de desenho de ${fmt(ST.deff,1)}, calculou-se a amostra de cada estrato separadamente, com correcção para população finita (Cochran, 1977)`}${ST.min>0?`, com um mínimo de ${fmt(ST.min)} entrevistas por estrato`:''}. A amostra efectiva é de ${fmt(R.net)} entrevistas; considerando uma taxa de resposta de ${fmt(ST.resp)}%, planeiam-se ${fmt(R.field)} entrevistas. A margem de erro esperada para o conjunto é de ±${fmt(R.overall*100,1)} pontos percentuais${H>1?`, variando entre ±${fmt(Math.min(...R.rows.map(r=>r.marg))*100,1)} e ±${fmt(maxM*100,1)} nos estratos`:''}.`;
  const en=`The sample was stratified by ${H} ${lvlEn}${scopeEn}, using population figures from the 2017 Census (INE). ${ST.target==='global'?`To estimate a prevalence of ${fmtEn(ST.p)}% with a margin of error of ±${fmtEn(ST.d,1)} percentage points at ${fmtEn(ST.conf*100)}% confidence and a design effect of ${fmtEn(ST.deff,1)}, a total sample of ${fmtEn(Math.ceil(R.nTot))} interviews was calculated (Cochran, 1977) and allocated ${ALLOC_EN[ST.alloc]}`:`So that each stratum achieves a margin of error of ±${fmtEn(ST.d,1)} percentage points at ${fmtEn(ST.conf*100)}% confidence with a design effect of ${fmtEn(ST.deff,1)}, each stratum's sample was calculated separately with a finite population correction (Cochran, 1977)`}${ST.min>0?`, with a minimum of ${fmtEn(ST.min)} interviews per stratum`:''}. The effective sample is ${fmtEn(R.net)} interviews; allowing for a ${fmtEn(ST.resp)}% response rate, ${fmtEn(R.field)} interviews are planned. The expected overall margin of error is ±${fmtEn(R.overall*100,1)} percentage points${H>1?`, ranging from ±${fmtEn(Math.min(...R.rows.map(r=>r.marg))*100,1)} to ±${fmtEn(maxM*100,1)} across strata`:''}.`;
  const body=R.rows.map(r=>`<tr><td>${esc(r.n)}${ST.level==='dist'&&ST.prov==='*'?`<span class="tagc">${esc(r.p)}</span>`:''}${r.minHit?'<span class="tagc">mínimo</span>':''}</td>
    <td class="num">${fmt(r.pop)}</td><td class="num">${fmt(r.W*100,1)}%</td>
    <td><input type="number" min="1" max="99" step="1" value="${fmt(r.ph*100)}" data-ph="${esc(r.key)}" aria-label="Prevalência esperada em ${esc(r.n)}"></td>
    <td class="num">${fmt(r.ns)}</td><td class="num"><b>${fmt(r.nField)}</b></td>
    <td class="num${r.marg>0.10?' bad':''}">±${fmt(r.marg*100,1)}<span class="mbar"><i style="width:${Math.min(100,r.marg/0.2*100).toFixed(0)}%"></i></span></td></tr>`).join('');
  $('st-results').innerHTML=`
   <div class="kpis">
    <div class="kpi main"><b>${fmt(R.field)}</b><small>entrevistas a planear</small></div>
    <div class="kpi"><b>±${fmt(R.overall*100,1)}</b><small>p.p. de margem no total</small></div>
    <div class="kpi"><b>${fmt(R.cost)}</b><small>MT de custo de campo</small></div>
    <div class="kpi"><b>${fmt(R.days)}</b><small>dias com ${fmt(ST.team)} inquiridores</small></div>
   </div>
   ${bad?`<div class="warn">${bad} de ${H} estratos ficam com margem acima de ±10 p.p. Se quer resultados fiáveis em cada estrato, use "Margem em cada estrato", a alocação pela raiz quadrada ou um mínimo por estrato maior.</div>`:''}
   ${R.minHits?`<div class="note">${R.minHits} estrato(s) subiram até ao mínimo de ${fmt(ST.min)} entrevistas.</div>`:''}
   <div><div class="sec-t">Alocação por estrato · ${H} ${lvl}<span class="actions"><button type="button" class="btn ghost" id="st-csv">Descarregar CSV</button></span></div>
    <div class="st-table"><table><thead><tr><th>Estrato</th><th>Pop. 2017</th><th>Peso</th><th>p (%)</th><th>Efectiva</th><th>A planear</th><th>Margem ±p.p.</th></tr></thead>
    <tbody>${body}</tbody><tfoot><tr><td>Total</td><td class="num">${fmt(R.N)}</td><td class="num">100%</td><td></td><td class="num">${fmt(R.net)}</td><td class="num">${fmt(R.field)}</td><td class="num">±${fmt(R.overall*100,1)}</td></tr></tfoot></table></div>
    <p class="note">Margem por estrato = z·√[p(1 − p)·deff / n · (1 − n/N)]. Margem total do estimador estratificado = z·√[Σ W²·p(1 − p)·deff / n · (1 − n/N)]. Orçamento: ${fmt(R.field)} entrevistas × ${fmt(ST.cost)} MT; ${fmt(R.enumDays,0)} dias-inquiridor a ${fmt(ST.perday)} entrevistas por dia.</p></div>
   <div><div class="sec-t">Texto para o protocolo
     <span style="display:flex;gap:8px;align-items:center"><span class="seg" role="group" aria-label="Idioma"><button type="button" data-sl="pt" aria-pressed="${ST.lang==='pt'}">PT</button><button type="button" data-sl="en" aria-pressed="${ST.lang==='en'}">EN</button></span><button type="button" class="btn" id="st-copy">Copiar</button></span></div>
     <p class="proto" id="st-proto">${esc(ST.lang==='pt'?pt:en)}</p></div>
   <div><div class="sec-t">Obras que inspiraram este cálculo</div><ul class="refs">${['cochran','kish','lohr','neyman'].map(k=>`<li>${REFS[k].f}</li>`).join('')}</ul></div>`;
  $('st-results').querySelectorAll('[data-ph]').forEach(el=>el.addEventListener('change',()=>{const v=parseFloat(el.value);if(v>0&&v<100){ST.pH[el.dataset.ph]=v;}else delete ST.pH[el.dataset.ph];calcStrata();}));
  $('st-results').querySelectorAll('[data-sl]').forEach(b=>b.addEventListener('click',()=>{ST.lang=b.dataset.sl;calcStrata();}));
  $('st-copy').addEventListener('click',e=>copyText(e.currentTarget,()=>$('st-proto').textContent,$('st-proto')));
  $('st-csv').addEventListener('click',()=>{
    const head=['estrato','provincia','populacao_2017','peso','p_esperada','amostra_efectiva','a_planear','margem_pp'];
    const lines=[head.join(';')].concat(R.rows.map(r=>[r.n,r.p,r.pop,(r.W).toFixed(4),(r.ph*100).toFixed(0),r.ns,r.nField,(r.marg*100).toFixed(2)].map(csvCell).join(';')));
    download('amostralab-alocacao.csv',lines.join('\n'));
  });
}

/* ================= quick data search ================= */
function initSearch(){
  const q=$('q'); q.addEventListener('input',renderSearch); renderSearch();
  const nD=DATA.districts.length;
  $('q-facts').innerHTML=`<span><b>${fmt(DATA.national)}</b> habitantes em 2017</span><span><b>11</b> províncias</span><span><b>${nD}</b> distritos e cidades</span>`;
}
function renderSearch(){
  const f=norm($('q').value.trim());
  let list, title;
  const provHit=f&&DATA.provinces.find(p=>norm(p.n)===f||norm(p.n).startsWith(f));
  if(!f){list=DATA.districts.slice().sort((a,b)=>b.pop-a.pop).slice(0,12);title='Os 12 distritos e cidades mais populosos';}
  else{list=DATA.districts.filter(x=>norm(x.n).includes(f)||norm(x.p).includes(f)||(x.from&&norm(x.from).includes(f))).sort((a,b)=>b.pop-a.pop);title=`${list.length} resultado(s)`;}
  let summary='';
  if(provHit){const ds=DATA.districts.filter(x=>x.p===provHit.n);summary=`<tr><td colspan="6" style="background:var(--accent-soft)"><b>${esc(provHit.n)}</b>: ${fmt(provHit.pop)} habitantes (${fmt(provHit.pop/DATA.national*100,1)}% do país), ${ds.length} distritos e cidades. <button type="button" class="link" data-addprov="${esc(provHit.n)}">Estratificar esta província por distrito</button></td></tr>`;}
  const max=Math.max(...list.map(x=>x.pop),1);
  $('q-results').innerHTML=list.length||summary?`<table><thead><tr><th>Distrito / cidade</th><th>Província</th><th class="num">População 2017</th><th class="num">% da província</th><th class="num">% do país</th><th></th></tr></thead><tbody>${summary}${list.slice(0,40).map(x=>{const pp=x.pop/DATA.provPop[x.p];return `<tr><td>${esc(x.n)}<span class="sub">${x.city?'Cidade / município':'Distrito'}${x.from?` · criado a partir de ${esc(x.from)}`:''}</span></td><td>${esc(x.p)}</td><td class="num">${fmt(x.pop)}<span class="share"><i style="width:${(x.pop/max*100).toFixed(0)}%"></i></span></td><td class="num">${fmt(pp*100,1)}%</td><td class="num">${fmt(x.pop/DATA.national*100,2)}%</td><td><button type="button" class="link" data-add="${esc('D:'+x.p+'|'+x.n)}" data-prov="${esc(x.p)}">+ Estrato</button></td></tr>`;}).join('')}</tbody></table>`:'<p class="note" style="padding:12px">Nenhum lugar com esse nome. Tente sem acentos ou só parte do nome.</p>';
  $('q-results').insertAdjacentHTML('afterbegin',`<div class="note" style="padding:10px 12px;border-bottom:1px solid var(--line)">${title}${list.length>40?' · a mostrar os 40 maiores':''}</div>`);
  $('q-results').querySelectorAll('[data-add]').forEach(b=>b.addEventListener('click',()=>{
    if(ST.level!=='dist'){ST.level='dist';ST.prov=b.dataset.prov;ST.sel=new Set();}
    else if(ST.prov!=='*'&&ST.prov!==b.dataset.prov){ST.prov='*';}
    ST.sel.add(b.dataset.add); ST.filter=''; $('st-search').value=''; syncSeg(); renderUnits(); calcStrata();
    b.textContent='Adicionado'; document.getElementById('estratificacao').scrollIntoView({behavior:'smooth',block:'start'});
  }));
  $('q-results').querySelectorAll('[data-addprov]').forEach(b=>b.addEventListener('click',()=>{
    ST.level='dist';ST.prov=b.dataset.addprov;selectAllInView();ST.filter='';$('st-search').value='';syncSeg();renderUnits();calcStrata();
    document.getElementById('estratificacao').scrollIntoView({behavior:'smooth',block:'start'});
  }));
}

/* ================= reproducible randomisation ================= */
function mulberry32(a){return function(){a|=0;a=a+0x6D2B79F5|0;let t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296;};}
const RZ={lang:'pt'};
(function initRz(){
  const ex=[];['Angoche','Monapo','Mogovolas'].forEach(dn=>{for(let i=1;i<=8;i++) ex.push(`EP ${dn} ${String(i).padStart(2,'0')}; ${dn}`);});
  $('rz-list').value=ex.join('\n'); $('rz-seed').value=20261117;
  ['rz-list','rz-seed','rz-arms'].forEach(id=>$(id).addEventListener('input',runRz));
  runRz();
})();
function runRz(){
  const lines=$('rz-list').value.split('\n').map(s=>s.trim()).filter(Boolean);
  const seed=parseInt($('rz-seed').value,10), arms=+$('rz-arms').value;
  if(!lines.length){$('rz-results').innerHTML='<div class="warn">Cole pelo menos uma unidade.</div>';return;}
  if(isNaN(seed)){$('rz-results').innerHTML='<div class="warn">Indique uma semente (número inteiro).</div>';return;}
  const units=lines.map((l,i)=>{const parts=l.split(/[;\t,]/);return {i,name:parts[0].trim(),stratum:(parts[1]||'').trim()||'Único'};});
  const strata=[...new Set(units.map(u=>u.stratum))];
  const rnd=mulberry32(seed), labels=arms===2?['Tratamento','Controlo']:['Tratamento A','Tratamento B','Controlo'], cls=arms===2?['T','C']:['T','T2','C'];
  let offset=0;
  strata.forEach(s=>{
    const g=units.filter(u=>u.stratum===s);
    for(let i=g.length-1;i>0;i--){const j=Math.floor(rnd()*(i+1));[g[i],g[j]]=[g[j],g[i]];}
    g.forEach((u,k)=>{u.arm=(k+offset)%arms;}); offset=(offset+g.length)%arms; // carry remainder so odd strata don't all favour the same arm
  });
  const counts=strata.map(s=>{const g=units.filter(u=>u.stratum===s);return {s,c:labels.map((_,a)=>g.filter(u=>u.arm===a).length)};});
  const tot=labels.map((_,a)=>units.filter(u=>u.arm===a).length);
  const pt=`A atribuição ${arms===2?'ao tratamento e ao controlo':'aos dois braços de tratamento e ao controlo'} foi feita por sorteio${strata.length>1?` estratificado por ${strata.length} estratos (${strata.join(', ')})`:''}, com ${units.length} unidades e a semente ${seed}, usando o gerador pseudo-aleatório Mulberry32 no AmostraLab. Resultado: ${labels.map((l,a)=>`${tot[a]} em ${l.toLowerCase()}`).join(', ')}. A aleatorização dentro de estratos melhora o equilíbrio entre grupos (Bruhn e McKenzie, 2009; Glennerster e Takavarasha, 2013). Com a mesma lista, pela mesma ordem, e a mesma semente, o sorteio pode ser reproduzido.`;
  const en=`Assignment to ${arms===2?'treatment and control':'two treatment arms and control'} was done by lottery${strata.length>1?` stratified by ${strata.length} strata (${strata.join(', ')})`:''}, with ${units.length} units and seed ${seed}, using the Mulberry32 pseudo-random generator in AmostraLab. Result: ${labels.map((l,a)=>`${tot[a]} in ${l.toLowerCase()}`).join(', ')}. Randomising within strata improves balance between groups (Bruhn and McKenzie, 2009; Glennerster and Takavarasha, 2013). The draw can be reproduced with the same list, in the same order, and the same seed.`;
  $('rz-results').innerHTML=`
   <div class="kpis">${labels.map((l,a)=>`<div class="kpi${a===0?' main':''}"><b>${tot[a]}</b><small>${l}</small></div>`).join('')}<div class="kpi"><b>${strata.length}</b><small>estrato(s)</small></div></div>
   <div><div class="sec-t">Resultado do sorteio<span class="actions"><button type="button" class="btn ghost" id="rz-csv">Descarregar CSV</button></span></div>
   <div class="st-table"><table><thead><tr><th>Unidade</th><th>Estrato</th><th>Grupo</th></tr></thead><tbody>${units.map(u=>`<tr><td>${esc(u.name)}</td><td style="text-align:left">${esc(u.stratum)}</td><td class="arm-${cls[u.arm]}" style="text-align:left">${labels[u.arm]}</td></tr>`).join('')}</tbody></table></div></div>
   <div><div class="sec-t">Equilíbrio por estrato</div><div class="st-table" style="max-height:none"><table><thead><tr><th>Estrato</th>${labels.map(l=>`<th>${l}</th>`).join('')}</tr></thead><tbody>${counts.map(c=>`<tr><td>${esc(c.s)}</td>${c.c.map(x=>`<td class="num">${x}</td>`).join('')}</tr>`).join('')}</tbody></table></div></div>
   <div><div class="sec-t">Texto para o protocolo<span style="display:flex;gap:8px;align-items:center"><span class="seg" role="group" aria-label="Idioma"><button type="button" data-rl="pt" aria-pressed="${RZ.lang==='pt'}">PT</button><button type="button" data-rl="en" aria-pressed="${RZ.lang==='en'}">EN</button></span><button type="button" class="btn" id="rz-copy">Copiar</button></span></div><p class="proto" id="rz-proto">${esc(RZ.lang==='pt'?pt:en)}</p></div>`;
  $('rz-results').querySelectorAll('[data-rl]').forEach(b=>b.addEventListener('click',()=>{RZ.lang=b.dataset.rl;runRz();}));
  $('rz-copy').addEventListener('click',e=>copyText(e.currentTarget,()=>$('rz-proto').textContent,$('rz-proto')));
  $('rz-csv').addEventListener('click',()=>download(`amostralab-sorteio-semente-${seed}.csv`,['unidade;estrato;grupo'].concat(units.map(u=>[u.name,u.stratum,labels[u.arm]].map(csvCell).join(';'))).join('\n')));
}

})();
