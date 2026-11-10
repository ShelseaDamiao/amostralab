/* AmostraLab — assistente de amostragem probabilística */
(function(){
const A=window.AL; if(!A) return;
const {normInv,fmt,fmtEn,esc,REFS,download,copyText,mulberry32,norm}=A;
const $=id=>document.getElementById(id);

Object.assign(REFS,{
  groves:{c:'Groves et al., 2009',f:'Groves, R. M., Fowler, F. J., Couper, M. P., Lepkowski, J. M., Singer, E. e Tourangeau, R. (2009). <i>Survey Methodology</i> (2.ª ed.). Wiley.',u:'População-alvo, base de amostragem e erros de cobertura e não resposta.',m:'Assistente'},
  valliant:{c:'Valliant et al., 2018',f:'Valliant, R., Dever, J. A. e Kreuter, F. (2018). <i>Practical Tools for Designing and Weighting Survey Samples</i> (2.ª ed.). Springer.',u:'Desenho multietápico, PPS e ponderação.',m:'Assistente'},
  un2005:{c:'Nações Unidas, 2005',f:'United Nations Statistics Division (2005). <i>Designing Household Survey Samples: Practical Guidelines</i>. Studies in Methods, Series F, No. 98. United Nations.',u:'Amostras domiciliares multietápicas em países em desenvolvimento.',m:'Assistente'},
  dhs2012:{c:'ICF International, 2012',f:'ICF International (2012). <i>Demographic and Health Survey Sampling and Household Listing Manual</i>. MEASURE DHS.',u:'Selecção PPS de áreas de enumeração e listagem de agregados.',m:'Assistente'},
  strobe:{c:'von Elm et al., 2007',f:'von Elm, E., Altman, D. G., Egger, M., Pocock, S. J., Gøtzsche, P. C. e Vandenbroucke, J. P. (2007). The Strengthening the Reporting of Observational Studies in Epidemiology (STROBE) statement. <i>The Lancet</i>, 370(9596), 1453–1457.',u:'Lista de itens a reportar sobre participantes e tamanho da amostra.',m:'Assistente'},
  lumley:{c:'Lumley, 2010',f:'Lumley, T. (2010). <i>Complex Surveys: A Guide to Analysis Using R</i>. Wiley.',u:'Análise de amostras complexas com pesos.',m:'Assistente'}
});
$('libtable').innerHTML=Object.values(REFS).map(r=>`<tr><td>${r.f}</td><td>${r.u}</td><td>${r.m}</td></tr>`).join('');

/* ---------------- designs ---------------- */
const D={
 aas:{n:'Aleatória simples',en:'simple random sampling',k:'AAS',
  def:'Cada unidade da base de amostragem tem a mesma probabilidade de ser seleccionada, e todas as amostras possíveis de tamanho n são igualmente prováveis.',
  defEn:'every unit in the sampling frame has the same probability of selection, and every possible sample of size n is equally likely',
  when:'Há uma lista completa e actualizada, e a população não está muito dispersa no terreno.',
  req:'Lista completa com todas as unidades numeradas de 1 a N.',
  pro:['É a base teórica de todos os outros desenhos','Análise simples, sem pesos'],
  con:['Exige lista completa','Em populações dispersas, o trabalho de campo fica caro','Subgrupos pequenos podem ficar com poucas unidades'],
  deff:'deff = 1 (referência)',refs:['cochran','lohr']},
 sist:{n:'Sistemática',en:'systematic sampling',k:'SIS',
  def:'Escolhe-se um ponto de partida aleatório e depois uma unidade em cada k da lista, com k = N / n (intervalo de amostragem).',
  defEn:'a random starting point is chosen and every k-th unit of the list is then selected, with k = N/n (sampling interval)',
  when:'Há uma lista ordenada ou um fluxo contínuo: registos, entradas sucessivas, casas ao longo de um percurso.',
  req:'Lista ordenada ou fluxo contínuo; N conhecido, pelo menos aproximadamente.',
  pro:['Muito fácil de aplicar no terreno','Se a lista estiver ordenada por uma variável relevante, funciona como estratificação implícita'],
  con:['Se a lista tiver um padrão que se repete de k em k, a amostra fica enviesada','Não há estimador exacto da variância; usa-se o da aleatória simples'],
  deff:'deff ≈ 1',refs:['cochran','lohr']},
 estrat:{n:'Estratificada',en:'stratified random sampling',k:'EST',
  def:'A população é dividida em grupos que não se sobrepõem (estratos), como províncias ou distritos, e faz-se uma selecção aleatória independente em cada um.',
  defEn:'the population is divided into non-overlapping groups (strata) and an independent random sample is drawn within each stratum',
  when:'Quer resultados fiáveis por subgrupo, ou os subgrupos são muito diferentes entre si.',
  req:'Saber a que estrato pertence cada unidade e o tamanho de cada estrato.',
  pro:['Garante que todos os estratos ficam representados','Mais precisa quando os estratos são homogéneos por dentro','Permite estimativas por estrato'],
  con:['Precisa de informação prévia para formar os estratos','Alocação não proporcional exige pesos na análise'],
  deff:'deff ≤ 1 com alocação proporcional',refs:['cochran','neyman','kish']},
 cong:{n:'Por conglomerados',en:'cluster sampling',k:'CON',
  def:'Sorteiam-se grupos naturais (conglomerados), como localidades, instituições ou áreas de enumeração, e inquirem-se as unidades dentro dos grupos sorteados.',
  defEn:'natural groups (clusters) such as villages, schools or enumeration areas are randomly selected and units within the selected clusters are surveyed',
  when:'Não há lista de unidades individuais mas há lista de grupos, ou a população está muito dispersa.',
  req:'Lista de conglomerados; listagem das unidades só nos conglomerados sorteados.',
  pro:['Reduz custos de deslocação','Só é preciso listar unidades nos conglomerados sorteados'],
  con:['Menos precisa: unidades do mesmo grupo parecem-se (ICC)','Precisa de uma amostra maior: deff = 1 + (m − 1)·ICC'],
  deff:'deff = 1 + (m − 1)·ICC',refs:['kish','lohr']},
 multi:{n:'Multietápica com PPS',en:'multi-stage sampling with probability proportional to size (PPS)',k:'MUL',
  def:'Selecção em etapas: primeiro áreas (unidades primárias) com probabilidade proporcional ao tamanho (PPS), depois um número fixo de unidades em cada área. Normalmente estratificada na primeira etapa.',
  defEn:'selection in stages: primary sampling units (areas) are drawn with probability proportional to size (PPS), followed by a fixed number of units within each selected area, usually with stratification at the first stage',
  when:'Inquéritos de grande escala, provinciais ou nacionais, sem lista de unidades individuais. É o desenho padrão dos grandes inquéritos por amostragem.',
  req:'Lista de áreas (por exemplo áreas de enumeração) com o seu tamanho; listagem das unidades nas áreas sorteadas.',
  pro:['Viável em grande escala sem lista de unidades individuais','PPS com número fixo de unidades por área torna a amostra autoponderada','Combina estratificação e conglomeração'],
  con:['A análise exige pesos e métodos para amostras complexas','deff típico entre 1,5 e 2,5','Mais etapas para documentar'],
  deff:'deff = 1 + (m − 1)·ICC, compensado em parte pela estratificação',refs:['kish','un2005','dhs2012','valliant']}
};
const UNIT={pessoa:['a pessoa','pessoas','person','people'],agregado:['o agregado familiar','agregados familiares','household','households'],instituicao:['a instituição ou organização','instituições','institution','institutions'],outra:['a unidade definida no estudo','unidades','unit','units']};

/* ---------------- state ---------------- */
const EX={step:1,title:'',pop:'',unit:'pessoa',incl:'',excl:'',
 nmode:'grande',narea:'N',npct:100,nval:'',frame:'lista',obj:'estimar',out:'prop',geo:'local',sub:'nao',design:'aas',
 p:50,src:'',sigma:10,src2:'',d:5,dMean:1,conf:0.95,resp:90,m:15,icc:0.05,stratify:false,
 cost:0,perday:6,team:10,seed:20260101,lang:'pt',c:''};
const W=JSON.parse(JSON.stringify(EX));
const ST={level:'prov',prov:'*',sel:new Set(),pH:{},target:'global',alloc:'prop',min:0,filter:''};
let DATA=null;

/* ---------------- helpers ---------------- */
const lc=t=>{t=String(t).trim();return t?t.charAt(0).toLowerCase()+t.slice(1)+(/[.!?]$/.test(t)?'':'.'):t;};
const isClu=()=>W.design==='cong'||W.design==='multi';
const isStrat=()=>W.design==='estrat'||(W.design==='multi'&&W.stratify);
const zOf=()=>normInv(1-(1-W.conf)/2);
const U=(i)=>UNIT[W.unit][i];
function areaOptions(){
  if(!DATA) return '';
  let h=`<option value="N">Moçambique (${fmt(DATA.national)})</option>`;
  DATA.provinces.forEach(p=>{h+=`<optgroup label="${esc(p.n)}"><option value="P:${esc(p.n)}">Província de ${esc(p.n)} (${fmt(p.pop)})</option>`+DATA.districts.filter(x=>x.p===p.n).map(x=>`<option value="D:${esc(p.n)}|${esc(x.n)}">${esc(x.n)} (${fmt(x.pop)})</option>`).join('')+'</optgroup>';});
  return h;
}
function areaPop(v){if(!DATA)return 0;if(v==='N')return DATA.national;if(v.startsWith('P:'))return DATA.provPop[v.slice(2)]||0;const [p,n]=v.slice(2).split('|');const x=DATA.districts.find(d=>d.p===p&&d.n===n);return x?x.pop:0;}
function areaName(v){if(v==='N')return 'Moçambique';if(v.startsWith('P:'))return 'província de '+v.slice(2);const [p,n]=v.slice(2).split('|');return `distrito de ${n} (${p})`;}
function getN(){
  if(isStrat()&&W.nmode==='censo') { const s=strataUnits(); if(s.length) return Math.round(s.reduce((a,u)=>a+u.pop,0)*elig()); }
  if(W.nmode==='censo') return Math.round(areaPop(W.narea)*elig());
  if(W.nmode==='manual') return parseInt(W.nval,10)||0;
  return 0;
}
const elig=()=>W.nmode==='censo'?(W.npct||100)/100:1;

/* ---------------- strata ---------------- */
const keyOf=u=>u.kind+':'+u.p+'|'+u.n;
function unitsInView(){
  if(!DATA) return [];
  if(ST.level==='prov') return DATA.provinces.map(p=>({kind:'P',p:p.n,n:p.n,pop:p.pop}));
  return DATA.districts.filter(x=>ST.prov==='*'||x.p===ST.prov).map(x=>({kind:'D',p:x.p,n:x.n,pop:x.pop}));
}
const strataUnits=()=>unitsInView().filter(u=>ST.sel.has(keyOf(u)));
const selectAllInView=()=>{ST.sel=new Set(unitsInView().map(keyOf));};
const ALLOC_H={prop:'Entrevistas na proporção da população. Óptima para o total, fraca para estratos pequenos.',raiz:'Proporcional à raiz quadrada da população (Kish). Reforça os estratos pequenos com pouca perda no total.',igual:'O mesmo número em cada estrato. Boa para comparar estratos.',neyman:'Proporcional à população × desvio-padrão (Neyman). Edite a prevalência de cada estrato na tabela do passo 5.'};
const TARGET_H={global:'Calcula a amostra para a margem pedida no total e distribui-a pelos estratos.',estrato:'Cada estrato recebe a amostra necessária para ter, sozinho, a margem pedida. A amostra total cresce muito.'};
function initStrataUI(){
  $('st-prov').innerHTML=`<option value="*">Todas as províncias</option>`+DATA.provinces.map(p=>`<option value="${esc(p.n)}">${esc(p.n)}</option>`).join('');
  $('st-prov').addEventListener('change',e=>{ST.prov=e.target.value;selectAllInView();renderUnits();update();});
  const seg=(id,k,after)=>$(id).querySelectorAll('button').forEach(b=>b.addEventListener('click',()=>{ST[k]=b.dataset.v;if(after)after();syncStrata();update();}));
  seg('st-level','level',()=>{selectAllInView();renderUnits();}); seg('st-target','target'); seg('st-alloc','alloc');
  $('st-search').addEventListener('input',e=>{ST.filter=e.target.value;renderUnits();});
  $('st-all').addEventListener('click',()=>{filtered().forEach(u=>ST.sel.add(keyOf(u)));renderUnits();update();});
  $('st-none').addEventListener('click',()=>{filtered().forEach(u=>ST.sel.delete(keyOf(u)));renderUnits();update();});
  $('st-min').value=ST.min; $('st-min').addEventListener('input',e=>{ST.min=parseFloat(e.target.value)||0;update();});
  selectAllInView(); syncStrata(); renderUnits();
}
function syncStrata(){
  [['st-level','level'],['st-target','target'],['st-alloc','alloc']].forEach(([id,k])=>$(id).querySelectorAll('button').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.v===ST[k]))));
  $('st-prov-wrap').hidden=ST.level!=='dist'; $('st-prov').value=ST.prov;
  $('st-alloc-h').textContent=ALLOC_H[ST.alloc]; $('st-target-h').textContent=TARGET_H[ST.target];
  $('st-alloc').style.opacity=ST.target==='estrato'?.45:1;
}
const filtered=()=>{const f=norm(ST.filter.trim());return unitsInView().filter(u=>!f||norm(u.n).includes(f)||norm(u.p).includes(f));};
function renderUnits(){
  const us=filtered(), all=unitsInView();
  $('st-units').innerHTML=us.length?us.map(u=>{const k=keyOf(u);return `<label><span><input type="checkbox" data-key="${esc(k)}"${ST.sel.has(k)?' checked':''}>${esc(u.n)}${ST.level==='dist'&&ST.prov==='*'?`<span class="tagc">${esc(u.p)}</span>`:''}</span><span class="pp">${fmt(u.pop)}</span></label>`;}).join(''):'<p class="note" style="padding:8px 10px;margin:0">Nenhuma unidade com esse nome.</p>';
  $('st-units').querySelectorAll('input').forEach(c=>c.addEventListener('change',()=>{c.checked?ST.sel.add(c.dataset.key):ST.sel.delete(c.dataset.key);count();update();}));
  count();
  function count(){$('st-count').textContent=`${all.filter(u=>ST.sel.has(keyOf(u))).length} de ${all.length} seleccionados`;}
}

/* ---------------- core calculation ---------------- */
function calc(){
  const z=zOf(), N=getN(), resp=(W.resp||100)/100, prop=W.out==='prop';
  const d=prop?W.d/100:W.dMean, p=W.p/100;
  if(!(d>0)) return {error:'Indique uma margem de erro maior que zero.'};
  if(prop&&!(p>0&&p<1)) return {error:'A prevalência esperada tem de estar entre 1 e 99%.'};
  if(!prop&&!(W.sigma>0)) return {error:'Indique o desvio-padrão esperado.'};
  const n0=prop?z*z*p*(1-p)/(d*d):z*z*W.sigma*W.sigma/(d*d);
  const nF=N>0?n0/(1+(n0-1)/N):n0;
  const deff=isClu()?1+(Math.max(1,W.m)-1)*W.icc:1;
  const nD=nF*deff;
  let R={z,N,n0,nF,deff,nD,resp,prop,d,p};
  if(isStrat()){
    const units=strataUnits();
    if(!units.length) return {...R,error:'Seleccione pelo menos um estrato no passo 4.'};
    const Nt=units.reduce((s,u)=>s+u.pop*elig(),0);
    const rows=units.map(u=>{const k=keyOf(u);const ph=(ST.pH[k]??W.p)/100;const Nh=u.pop*elig();return {...u,key:k,ph,Nh,W:Nh/Nt,S:prop?Math.sqrt(ph*(1-ph)):1};});
    if(ST.target==='global'||!prop){
      const w=r=>ST.alloc==='prop'?r.Nh:ST.alloc==='raiz'?Math.sqrt(r.Nh):ST.alloc==='igual'?1:r.Nh*r.S;
      const sw=rows.reduce((s,r)=>s+w(r),0); rows.forEach(r=>{r.nRaw=nD*w(r)/sw;});
    } else rows.forEach(r=>{const a=z*z*r.ph*(1-r.ph)/(d*d);r.nRaw=a/(1+(a-1)/r.Nh)*deff;});
    rows.forEach(r=>{r.ns=Math.min(Math.floor(r.Nh),Math.max(ST.min,1,Math.ceil(r.nRaw)));r.minHit=r.ns>Math.ceil(r.nRaw);
      r.nField=Math.ceil(r.ns/resp); if(isClu()){r.cl=Math.max(2,Math.ceil(r.nField/W.m));r.minCl=r.cl*W.m>Math.ceil(r.ns/resp)&&Math.ceil(Math.ceil(r.ns/resp)/W.m)<2;r.nField=r.cl*W.m;}
      r.marg=prop?z*Math.sqrt(r.ph*(1-r.ph)*deff/r.ns*Math.max(0,1-r.ns/r.Nh)):null;});
    const V=prop?rows.reduce((s,r)=>s+r.W*r.W*r.ph*(1-r.ph)*deff/r.ns*Math.max(0,1-r.ns/r.Nh),0):0;
    R={...R,rows,net:rows.reduce((s,r)=>s+r.ns,0),field:rows.reduce((s,r)=>s+r.nField,0),overall:prop?z*Math.sqrt(V):null,clusters:isClu()?rows.reduce((s,r)=>s+r.cl,0):null,N:Math.round(Nt)};
  } else {
    let field=Math.ceil(nD/resp), clusters=null;
    if(isClu()){clusters=Math.max(2,Math.ceil(field/W.m));field=clusters*W.m;}
    R={...R,net:Math.ceil(nD),field,clusters,overall:d};
  }
  R.cost=R.field*(W.cost||0); R.enumDays=R.field/(W.perday||1); R.days=Math.ceil(R.enumDays/(W.team||1));
  if(R.N>0&&R.field>R.N) R.warn='A amostra é maior do que a população elegível. Considere um censo (inquirir todos).';
  return R;
}

/* ---------------- recommendation ---------------- */
function recommend(){
  if(W.frame==='nenhuma') return {d:'multi',why:'Sem base de amostragem, o caminho habitual é usar uma lista de áreas (por exemplo, as áreas de enumeração do censo) como base da primeira etapa e listar as unidades nas áreas sorteadas. Isso conduz a um desenho multietápico.'};
  if(W.frame==='areas') return (W.geo==='local'&&W.sub==='nao')?{d:'cong',why:'Tem lista de grupos, não de pessoas, e a área é pequena. Sortear conglomerados e inquirir dentro deles é o mais simples.'}:{d:'multi',why:'Tem lista de áreas, a população está dispersa'+(W.sub==='sim'?' e precisa de resultados por subgrupo':'')+'. O desenho multietápico estratificado com PPS é o padrão internacional para este caso (Nações Unidas, 2005).'};
  if(W.sub==='sim') return {d:'estrat',why:'Tem lista completa e precisa de resultados fiáveis por subgrupo. A estratificação garante amostra suficiente em cada um.'};
  if(W.geo==='ampla') return {d:'estrat',why:'Tem lista completa mas a população está muito dispersa. Estratificar por área geográfica melhora a precisão e organiza o trabalho de campo.'};
  return {d:'aas',why:'Tem lista completa e a população é compacta. A aleatória simples é o desenho de referência; a sistemática é uma alternativa prática se a lista estiver ordenada.'};
}

/* ---------------- illustration ---------------- */
function dots(key){
  const cols=16,rows=8,r=mulberry32(7),cell=14,W0=cols*cell,H0=rows*cell; let s='';
  const pick=new Set();
  if(key==='aas'){while(pick.size<16)pick.add(Math.floor(r()*cols*rows));}
  if(key==='sist'){for(let i=3;i<cols*rows;i+=8)pick.add(i);}
  if(key==='estrat'){for(let b=0;b<4;b++){const got=new Set();while(got.size<4)got.add(b*32+Math.floor(r()*32));got.forEach(x=>pick.add(x));}for(let b=0;b<4;b++)s+=`<rect x="0" y="${b*2*cell}" width="${W0}" height="${2*cell}" class="${b%2?'band2':'band'}"/>`;}
  const clusters=[];
  if(key==='cong'||key==='multi'){for(let cy=0;cy<2;cy++)for(let cx=0;cx<4;cx++)clusters.push({cx,cy});
    const sel=key==='cong'?[1,4,6]:[0,2,5,7];
    sel.forEach(ci=>{const c=clusters[ci];const cells=[];for(let y=0;y<4;y++)for(let x=0;x<4;x++)cells.push((c.cy*4+y)*cols+c.cx*4+x);
      if(key==='cong')cells.forEach(x=>pick.add(x));else{const g=new Set();while(g.size<4)g.add(cells[Math.floor(r()*16)]);g.forEach(x=>pick.add(x));}
      s+=`<rect x="${c.cx*4*cell+1}" y="${c.cy*4*cell+1}" width="${4*cell-2}" height="${4*cell-2}" rx="4" class="clu"/>`;});
    if(key==='multi'){s+=`<line x1="0" y1="${4*cell}" x2="${W0}" y2="${4*cell}" class="sline"/>`;}
    clusters.forEach(c=>{s+=`<rect x="${c.cx*4*cell+1}" y="${c.cy*4*cell+1}" width="${4*cell-2}" height="${4*cell-2}" rx="4" class="cluo"/>`;});
  }
  for(let i=0;i<cols*rows;i++){const x=(i%cols)*cell+cell/2,y=Math.floor(i/cols)*cell+cell/2;s+=`<circle cx="${x}" cy="${y}" r="${pick.has(i)?4.2:2.6}" class="${pick.has(i)?'on':'off'}"/>`;}
  return `<svg viewBox="0 0 ${W0} ${H0}" class="dots" role="img" aria-label="Ilustração: ${esc(D[key].n)}">${s}</svg>`;
}

/* ---------------- rendering ---------------- */
function bindText(id,k,cb){const el=$(id);el.value=W[k];el.addEventListener('input',()=>{W[k]=el.value;cb&&cb();renderSummary();});}
function bindNum(id,k,rangeId){const el=$(id),rg=rangeId&&$(rangeId);el.value=W[k];if(rg)rg.value=W[k];
  const set=v=>{W[k]=parseFloat(v);update();};
  el.addEventListener('input',()=>{if(rg)rg.value=el.value;set(el.value);});
  if(rg)rg.addEventListener('input',()=>{el.value=rg.value;set(rg.value);});}
function bindChoice(id,k,after){const box=$(id);box.querySelectorAll('button').forEach(b=>b.addEventListener('click',()=>{W[k]=b.dataset.v;syncChoices();after&&after();update();}));}
function syncChoices(){
  [['w-frame','frame'],['w-obj','obj'],['w-out','out'],['w-geo','geo'],['w-sub','sub'],['w-nmode','nmode']].forEach(([id,k])=>$(id).querySelectorAll('button').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.v===W[k]))));
  $('w-n-censo').hidden=W.nmode!=='censo'; $('w-n-manual').hidden=W.nmode!=='manual';
  $('w-obj-warn').hidden=W.obj!=='comparar';
}
function renderDesigns(){
  const rec=recommend();
  $('w-rec').innerHTML=`<div class="recbox"><span class="kicker">Recomendado para o seu caso</span><b>${D[rec.d].n}</b><p>${rec.why}</p>${W.design!==rec.d?`<button type="button" class="btn ghost" id="w-use-rec">Usar a recomendação</button>`:'<span class="note">Está a usar o desenho recomendado.</span>'}</div>`;
  const ub=$('w-use-rec'); if(ub) ub.addEventListener('click',()=>{W.design=rec.d;if(rec.d==='multi')W.stratify=true;update();renderDesigns();});
  $('w-designs').innerHTML=Object.entries(D).map(([k,x])=>`<article class="design${W.design===k?' sel':''}${rec.d===k?' rec':''}">
    <div class="dhead">${dots(k)}<div><h4>${x.n}${rec.d===k?' <span class="pill">recomendado</span>':''}</h4><p>${x.def}</p></div></div>
    <dl class="dmeta"><div><dt>Quando usar</dt><dd>${x.when}</dd></div><div><dt>O que precisa</dt><dd>${x.req}</dd></div>
    <div><dt>Vantagens</dt><dd><ul>${x.pro.map(t=>`<li>${t}</li>`).join('')}</ul></dd></div><div><dt>Limitações</dt><dd><ul>${x.con.map(t=>`<li>${t}</li>`).join('')}</ul></dd></div></dl>
    <div class="dfoot"><span class="mono">${x.deff}</span><button type="button" class="btn${W.design===k?'':' ghost'}" data-design="${k}">${W.design===k?'Seleccionado':'Escolher este desenho'}</button></div></article>`).join('');
  $('w-designs').querySelectorAll('[data-design]').forEach(b=>b.addEventListener('click',()=>{W.design=b.dataset.design;update();renderDesigns();}));
}
function renderParams(R){
  const prop=W.out==='prop';
  $('w-par-prop').hidden=!prop; $('w-par-mean').hidden=prop; $('w-d-u').textContent=prop?'p.p.':'unidades';
  $('w-d').value=prop?W.d:W.dMean; $('w-d-r').hidden=!prop; if(prop)$('w-d-r').value=W.d;
  $('w-par-clu').hidden=!isClu(); $('w-strat-wrap').hidden=W.design!=='multi'; $('w-stratify').checked=W.stratify;
  $('w-strata').hidden=!isStrat();
  $('w-par-lead').textContent=`Desenho escolhido: ${D[W.design].n}${W.design==='multi'&&W.stratify?' (estratificada)':''}. Cada campo explica o que significa e que valor usar quando não sabe.`;
  const cf=fmt(W.conf*100);
  $('w-d-explain').innerHTML=prop?`<b>Em linguagem simples:</b> se a amostra der ${fmt(W.p)}%, o valor real na população estará entre <b>${fmt(Math.max(0,W.p-W.d),1)}% e ${fmt(Math.min(100,W.p+W.d),1)}%</b>, com ${cf}% de confiança. Ou seja, em ${cf} de cada 100 amostras possíveis, o intervalo assim construído conteria o valor real.`
    :`<b>Em linguagem simples:</b> a média estimada ficará a menos de ${fmt(W.dMean,2)} unidades do valor real, com ${cf}% de confiança.`;
  if(isClu()) $('w-deff-out').innerHTML=`<b>Efeito de desenho:</b> deff = 1 + (${fmt(W.m)} − 1) × ${fmt(W.icc,2)} = <b>${fmt(1+(W.m-1)*W.icc,2)}</b>. A amostra tem de ser ${fmt(1+(W.m-1)*W.icc,2)} vezes maior do que numa aleatória simples. Menos unidades por conglomerado e mais conglomerados reduzem o deff.`;
  $('w-preview').innerHTML=R.error?`<div class="warn">${R.error}</div>`:`<div class="pv"><span class="pv-n">${fmt(R.field)}</span><span>${U(1)} a seleccionar${R.clusters?` em <b>${fmt(R.clusters)}</b> conglomerados`:''}${R.rows?` · ${R.rows.length} estratos`:''}</span><button type="button" class="btn" data-go="5">Ver detalhe</button></div>`;
  const g=$('w-preview').querySelector('[data-go]'); if(g) g.addEventListener('click',()=>go(5));
}
function stepsList(R){
  const L=[{l:'Amostra base, n₀',n:R.prop?`z²·p(1 − p)/d² com z = ${fmt(R.z,3)}, p = ${fmt(W.p)}%, d = ${fmt(W.d,1)} p.p.`:`z²·σ²/d² com σ = ${fmt(W.sigma,2)}, d = ${fmt(W.dMean,2)}`,v:R.n0},
   {l:'Correcção para população finita',n:R.N>0?`N = ${fmt(R.N)} ${U(1)} elegíveis`:'Não aplicada (N desconhecido ou muito grande)',v:R.nF}];
  if(isClu()) L.push({l:'× Efeito de desenho',n:`1 + (m − 1)·ICC = ${fmt(R.deff,2)}`,v:R.nD});
  L.push({l:'÷ Taxa de resposta'+(isClu()?' e arredondamento a conglomerados completos':''),n:`${fmt(W.resp)}% de resposta${isClu()?`; ${fmt(W.m)} por conglomerado`:''}${R.rows?`; soma dos ${R.rows.length} estratos`:''}`,v:R.field,final:true});
  const mx=Math.max(...L.map(s=>s.v));
  return `<div class="steps">${L.map(s=>`<div class="step${s.final?' final':''}"><span class="sl">${s.l}</span><span class="sv">${fmt(s.v,s.final?0:1)}</span><span class="sn">${s.n}</span><span class="bar"><i style="width:${(s.v/mx*100).toFixed(1)}%"></i></span></div>`).join('')}</div>`;
}
function renderResult(R){
  if(R.error){$('w-result').innerHTML=`<div class="warn">${R.error}</div>`;$('w-select').innerHTML='';return;}
  const prop=R.prop;
  let strata='';
  if(R.rows){
    strata=`<div class="sec-t">Alocação por estrato<span class="actions"><button type="button" class="btn ghost" id="w-st-csv">Descarregar CSV</button></span></div>
    <div class="st-table"><table><thead><tr><th>Estrato</th><th>N elegível</th><th>Peso</th>${prop?'<th>p (%)</th>':''}<th>Efectiva</th><th>A seleccionar</th>${isClu()?'<th>Conglom.</th>':''}${prop?'<th>Margem ±p.p.</th>':''}</tr></thead><tbody>
    ${R.rows.map(r=>`<tr><td>${esc(r.n)}${ST.level==='dist'&&ST.prov==='*'?`<span class="tagc">${esc(r.p)}</span>`:''}${r.minHit?'<span class="tagc">mínimo</span>':''}</td><td class="num">${fmt(r.Nh)}</td><td class="num">${fmt(r.W*100,1)}%</td>
      ${prop?`<td><input type="number" min="1" max="99" step="1" value="${fmt(r.ph*100)}" data-ph="${esc(r.key)}" aria-label="Prevalência em ${esc(r.n)}"></td>`:''}<td class="num">${fmt(r.ns)}</td><td class="num"><b>${fmt(r.nField)}</b></td>${isClu()?`<td class="num">${fmt(r.cl)}</td>`:''}
      ${prop?`<td class="num${r.marg>0.10?' bad':''}">±${fmt(r.marg*100,1)}<span class="mbar"><i style="width:${Math.min(100,r.marg/0.2*100).toFixed(0)}%"></i></span></td>`:''}</tr>`).join('')}
    </tbody><tfoot><tr><td>Total</td><td class="num">${fmt(R.N)}</td><td class="num">100%</td>${prop?'<td></td>':''}<td class="num">${fmt(R.net)}</td><td class="num">${fmt(R.field)}</td>${isClu()?`<td class="num">${fmt(R.clusters)}</td>`:''}${prop?`<td class="num">±${fmt(R.overall*100,1)}</td>`:''}</tr></tfoot></table></div>
    ${isClu()&&R.rows.some(r=>r.minCl)?`<div class="note">Cada estrato recebe pelo menos 2 conglomerados, o mínimo necessário para estimar a variância (Kish, 1965).</div>`:''}${prop&&W.sub==='sim'&&R.rows.some(r=>r.marg>0.10)?`<div class="warn">${R.rows.filter(r=>r.marg>0.10).length} estrato(s) ficam com margem acima de ±10 p.p. Para estimativas fiáveis em cada estrato, escolha "Margem em cada estrato" no passo 4, a alocação pela raiz quadrada ou aumente o mínimo por estrato.</div>`:''}`;
  }
  $('w-result').innerHTML=`
   <div class="kpis" style="margin-bottom:18px">
    <div class="kpi main"><b>${fmt(R.field)}</b><small>${U(1)} a seleccionar</small></div>
    <div class="kpi"><b>${R.clusters?fmt(R.clusters):fmt(R.net)}</b><small>${R.clusters?'conglomerados (áreas) a sortear':U(1)+' que se espera que respondam'}</small></div>
    <div class="kpi"><b>±${prop?fmt((R.rows?R.overall:R.d)*100,1)+' p.p.':fmt(R.d,2)}</b><small>margem de erro no total</small></div>
    <div class="kpi"><b>${fmt(R.days)}</b><small>dias de campo · ${W.cost>0?fmt(R.cost)+' MT':'custo por definir'}</small></div>
   </div>
   ${R.warn?`<div class="warn">${R.warn}</div>`:''}
   <div><div class="sec-t">Do cálculo ao terreno</div>${stepsList(R)}</div>
   ${strata}${R.rows&&prop&&W.sub!=='sim'?'<p class="note">Como o objectivo é o total, a estratificação serve para melhorar a precisão e garantir cobertura de todos os distritos. As margens por estrato são largas e não devem ser usadas para conclusões por distrito.</p>':''}`;
  $('w-result').querySelectorAll('[data-ph]').forEach(el=>el.addEventListener('change',()=>{const v=parseFloat(el.value);if(v>0&&v<100)ST.pH[el.dataset.ph]=v;else delete ST.pH[el.dataset.ph];update();}));
  const csv=$('w-st-csv'); if(csv) csv.addEventListener('click',()=>download('amostralab-alocacao.csv',['estrato;provincia;N_elegivel;peso;p;efectiva;a_seleccionar;conglomerados;margem_pp'].concat(R.rows.map(r=>[r.n,r.p,Math.round(r.Nh),r.W.toFixed(4),(r.ph*100).toFixed(0),r.ns,r.nField,r.cl??'',r.marg!=null?(r.marg*100).toFixed(2):''].join(';'))).join('\n')));
  renderSelection(R);
}

/* ---------------- selection procedures ---------------- */
function floydSample(N,n,rnd){const s=new Set();for(let j=N-n+1;j<=N;j++){const t=1+Math.floor(rnd()*j);s.add(s.has(t)?j:t);}return [...s].sort((a,b)=>a-b);}
function exampleEAs(){const r=mulberry32(11);const out=[];for(let i=1;i<=40;i++)out.push(`AE ${String(i).padStart(3,'0')}; ${80+Math.floor(r()*170)}`);return out.join('\n');}
let eaText=null;
function renderSelection(R){
  const seedField=`<div class="field" style="max-width:220px"><label for="w-seed">Semente do sorteio</label><div class="ctl"><input type="number" id="w-seed" value="${W.seed}" step="1"></div></div>`;
  const weightBox=(f,txt)=>`<div class="explain"><b>Peso amostral:</b> <span class="mono">${f}</span>. ${txt}</div>`;
  let h='<div class="sec-t" style="margin-top:20px">Como seleccionar a amostra</div>';
  const N=getN(), n=R.field;
  if(W.design==='aas'||W.design==='sist'){
    if(!(N>0)){$('w-select').innerHTML=h+`<div class="warn">Para sortear é preciso saber N, o número de unidades na lista. Indique-o no passo 1.</div>`;return;}
    if(n>=N){$('w-select').innerHTML=h+`<div class="warn">A amostra cobre toda a população: inquira todas as unidades.</div>`;return;}
    const rnd=mulberry32(W.seed|0);
    if(W.design==='aas'){
      const s=floydSample(N,n,rnd);
      h+=`<ol class="proc"><li>Numere as ${fmt(N)} unidades da lista de 1 a ${fmt(N)}, pela ordem em que aparecem.</li><li>Sorteie ${fmt(n)} números diferentes entre 1 e ${fmt(N)} (sem reposição). A lista abaixo foi gerada com a semente indicada e pode ser reproduzida.</li><li>Inquira as unidades com esses números. Não substitua quem não responde: a taxa de resposta já foi considerada.</li></ol>
      ${seedField}<div class="numgrid">${s.slice(0,80).map(x=>`<span>${x}</span>`).join('')}${s.length>80?`<span class="more">+ ${fmt(s.length-80)}</span>`:''}</div>
      <div class="actions"><button type="button" class="btn ghost" id="w-sel-csv">Descarregar os ${fmt(n)} números</button></div>
      ${weightBox(`w = N / n = ${fmt(N)} / ${fmt(R.net)} = ${fmt(N/R.net,2)}`,'Todas as unidades têm o mesmo peso: a amostra é autoponderada.')}`;
      $('w-select').innerHTML=h; $('w-sel-csv').addEventListener('click',()=>download(`amostralab-aas-semente-${W.seed}.csv`,['ordem;numero_seleccionado'].concat(s.map((x,i)=>`${i+1};${x}`)).join('\n')));
    } else {
      const k=N/n, r0=rnd()*k, s=[];for(let i=0;i<n;i++)s.push(Math.floor(r0+i*k)+1);
      h+=`<ol class="proc"><li>Ordene a lista (por exemplo por bairro ou por idade, o que dá uma estratificação implícita) e numere de 1 a ${fmt(N)}.</li><li>Intervalo de amostragem: k = N / n = ${fmt(N)} / ${fmt(n)} = <b>${fmt(k,3)}</b>.</li><li>Início aleatório entre 0 e k: <b>${fmt(r0,3)}</b> (semente ${W.seed}). Seleccione as posições ⌊início + i·k⌋ + 1, para i = 0, 1, …, n − 1.</li><li>Verifique que a lista não tem um padrão que se repete de ${fmt(Math.round(k))} em ${fmt(Math.round(k))} unidades.</li></ol>
      ${seedField}<div class="numgrid">${s.slice(0,60).map(x=>`<span>${x}</span>`).join('')}${s.length>60?`<span class="more">+ ${fmt(s.length-60)}</span>`:''}</div>
      <div class="actions"><button type="button" class="btn ghost" id="w-sel-csv">Descarregar as ${fmt(n)} posições</button></div>
      ${weightBox(`w = k = ${fmt(k,2)}`,'Igual para todas as unidades.')}`;
      $('w-select').innerHTML=h; $('w-sel-csv').addEventListener('click',()=>download(`amostralab-sistematica-semente-${W.seed}.csv`,['ordem;posicao'].concat(s.map((x,i)=>`${i+1};${x}`)).join('\n')));
    }
  } else if(W.design==='estrat'){
    h+=`<ol class="proc"><li>Separe a base de amostragem em ${R.rows.length} listas, uma por estrato.</li><li>Em cada estrato, sorteie o número de unidades indicado na coluna "A seleccionar", por aleatória simples ou sistemática, com sorteios independentes.</li><li>Guarde a semente usada em cada estrato para que o sorteio possa ser auditado.</li></ol>
    ${weightBox('w<sub>h</sub> = N<sub>h</sub> / n<sub>h</sub>',ST.alloc==='prop'?'Com alocação proporcional, os pesos são quase iguais.':'Como a alocação não é proporcional, os pesos variam entre estratos e têm de ser usados na análise (Lumley, 2010).')}`;
    $('w-select').innerHTML=h;
  } else {
    if(eaText===null) eaText=exampleEAs();
    const cDefault=R.rows?Math.max(...R.rows.map(r=>r.cl)):R.clusters;
    const c=Math.max(1,parseInt(W.c,10)||Math.min(cDefault,20));
    const lines=eaText.split('\n').map(s=>s.trim()).filter(Boolean).map(l=>{const p=l.split(/[;\t,]/);return {name:p[0].trim(),M:parseFloat(p[1])};}).filter(x=>x.name&&x.M>0);
    let tool='';
    if(lines.length<=c){tool=`<div class="warn">A lista tem ${lines.length} áreas e quer sortear ${c}. Cole uma lista maior ou reduza o número.</div>`;}
    else{
      const M=lines.reduce((s,x)=>s+x.M,0), I=M/c, rnd=mulberry32(W.seed|0), start=rnd()*I;
      let cum=0, j=0; const sel=[];
      lines.forEach(x=>{const lo=cum;cum+=x.M;while(j<c&&start+j*I<cum&&start+j*I>=lo){sel.push(x);j++;}});
      const mm=W.design==='multi'?W.m:null;
      tool=`<div class="st-table" style="max-height:320px"><table><thead><tr class="nocase"><th>Área sorteada</th><th>Tamanho (M<sub>i</sub>)</th><th>π₁ = c·M<sub>i</sub>/M</th>${mm?'<th>π₂ = m/M<sub>i</sub></th><th>Peso</th>':''}</tr></thead><tbody>${sel.map(x=>{const p1=Math.min(1,c*x.M/M),p2=mm?Math.min(1,mm/x.M):1;return `<tr><td>${esc(x.name)}${c*x.M>=M?'<span class="tagc">certeza</span>':''}</td><td class="num">${fmt(x.M)}</td><td class="num">${fmt(p1,4)}</td>${mm?`<td class="num">${fmt(p2,4)}</td><td class="num">${fmt(1/(p1*p2),1)}</td>`:''}</tr>`;}).join('')}</tbody></table></div>
      <p class="note">${sel.some(x=>c*x.M>=M)?'Há áreas maiores do que o intervalo I: devem ser incluídas com certeza (π₁ = 1), retiradas da lista, e o sorteio repetido com as restantes. ':''}Total M = ${fmt(M)} unidades · intervalo I = M / c = ${fmt(I,2)} · início aleatório = ${fmt(start,2)} (semente ${W.seed}).${mm?` Com PPS e ${mm} unidades fixas por área, o peso é igual em todas as áreas, M/(c·m) = ${fmt(M/(c*mm),1)}: a amostra é autoponderada dentro do estrato.`:''}</p>`;
    }
    h+=`<ol class="proc">
      <li><b>1.ª etapa: áreas.</b> ${R.rows?`Em cada um dos ${R.rows.length} estratos, sorteie o número de conglomerados da coluna "Conglom."`:`Sorteie ${fmt(R.clusters)} conglomerados`} ${W.design==='multi'?'com probabilidade proporcional ao tamanho (PPS sistemática), usando o tamanho de cada área.':'por aleatória simples a partir da lista de conglomerados.'}</li>
      <li><b>Listagem.</b> Nas áreas sorteadas, actualize a lista de todas as unidades antes da recolha de dados.</li>
      <li><b>2.ª etapa: ${U(1)}.</b> Em cada área, seleccione ${fmt(W.m)} ${U(1)} por amostragem sistemática a partir da listagem.</li>
      <li><b>Elegibilidade.</b> Aplique os critérios de inclusão e de exclusão a cada unidade seleccionada. Se uma unidade contiver mais de um elemento elegível, inclua todos ou sorteie um (e registe-o, porque altera os pesos).</li></ol>
      <div class="subbox"><div class="sec-t">Simulador de selecção PPS sistemática</div>
      <p class="note" style="margin-top:0">Cole a lista de áreas (ou outros conglomerados) de um estrato, com o tamanho de cada uma (nome ; tamanho). O exemplo usa áreas fictícias.</p>
      <div class="row2"><textarea id="w-ea" class="short" spellcheck="false">${esc(eaText)}</textarea>
      <div>${seedField}<div class="field" style="max-width:220px"><label for="w-c">Áreas a sortear (c)</label><div class="ctl"><input type="number" id="w-c" min="1" step="1" value="${c}"></div></div></div></div>${tool}</div>
      ${weightBox(W.design==='multi'?'w = 1 / (π₁ · π₂) = M / (c · m)':'w = C / c (todos os inquiridos no conglomerado)','Na análise, declare os estratos, os conglomerados e os pesos (por exemplo svyset no Stata ou o pacote survey no R; Lumley, 2010).')}`;
    $('w-select').innerHTML=h;
    $('w-ea').addEventListener('change',e=>{eaText=e.target.value;renderSelection(R);});
    $('w-c').addEventListener('change',e=>{W.c=e.target.value;renderSelection(R);});
  }
  const sd=$('w-seed'); if(sd) sd.addEventListener('change',e=>{W.seed=parseInt(e.target.value,10)||0;renderSelection(R);});
}

/* ---------------- protocol ---------------- */
function protocol(R,lang){
  const pt=lang==='pt', F=pt?fmt:fmtEn, x=D[W.design], prop=W.out==='prop', cf=F(W.conf*100);
  const strat=isStrat(), lvl=ST.level==='prov'?(pt?'províncias':'provinces'):(pt?'distritos':'districts');
  const sec=[];
  const H=(t)=>`<h5>${t}</h5>`;
  if(pt){
    const T=(v,ph)=>v&&v.trim()?esc(v.trim()):`<mark>[${ph}]</mark>`;
    sec.push(H('1. População do estudo')+`<p>${T(W.title,'objectivo do estudo')}</p><p>A população-alvo do estudo é constituída por: ${W.pop.trim()?esc(lc(W.pop)):'<mark>[definir a população-alvo]</mark>'} A unidade de análise é ${U(0)}.</p><p><b>Critérios de inclusão:</b> ${T(W.incl,'definir')}</p><p><b>Critérios de exclusão:</b> ${T(W.excl,'definir')}</p>${R.N>0?`<p>O tamanho da população elegível foi estimado em ${fmt(R.N)} ${U(1)}${W.nmode==='censo'&&!strat?`, aplicando ${fmt(W.npct,1)}% à população da ${areaName(W.narea)} no Censo 2017 (INE)`:W.nmode==='censo'?`, aplicando ${fmt(W.npct,1)}% à população dos estratos no Censo 2017 (INE)`:''}.</p>`:''}`);
    sec.push(H('2. Base de amostragem')+`<p>${W.frame==='lista'?'A base de amostragem é a lista completa das unidades da população-alvo, verificada quanto a duplicações e omissões antes do sorteio.':W.frame==='areas'?'A base de amostragem da primeira etapa é a lista de áreas (ou outros conglomerados) com o respectivo tamanho. Nas áreas seleccionadas, será feita uma listagem actualizada das unidades, que servirá de base para a segunda etapa (Groves et al., 2009).':'Não existe uma lista prévia das unidades. Será usada uma lista de áreas (por exemplo, as áreas de enumeração do INE) como base da primeira etapa, seguida da listagem completa das unidades nas áreas seleccionadas (Nações Unidas, 2005).'}</p>`);
    sec.push(H('3. Desenho amostral')+`<p>Será usada uma amostragem probabilística ${x.n.toLowerCase()}${strat&&W.design==='multi'?', estratificada por '+lvl:''}: ${x.def.charAt(0).toLowerCase()+x.def.slice(1)} ${strat?`Os estratos são ${R.rows.length} ${lvl}${ST.level==='dist'&&ST.prov!=='*'?' da província de '+ST.prov:''}.`:''} A escolha justifica-se porque ${x.when.charAt(0).toLowerCase()+x.when.slice(1)}</p>`);
    let calc=`<p>O tamanho da amostra foi calculado para estimar ${prop?`uma proporção esperada de ${fmt(W.p)}% com uma margem de erro de ±${fmt(W.d,1)} pontos percentuais`:`uma média com uma margem de erro de ±${fmt(W.dMean,2)} unidades, assumindo um desvio-padrão de ${fmt(W.sigma,2)}`} e ${cf}% de confiança (Cochran, 1977), o que dá uma amostra base de ${Math.ceil(R.n0)}. ${(prop?W.src:W.src2)?`O valor esperado baseia-se em: ${esc(prop?W.src:W.src2)}. `:''}${R.N>0?`Aplicou-se a correcção para população finita (N = ${fmt(R.N)}), obtendo-se ${Math.ceil(R.nF)}. `:''}${isClu()?`Para ter em conta a conglomeração, aplicou-se um efeito de desenho de ${fmt(R.deff,2)}, calculado como 1 + (m − 1)·ICC com ${fmt(W.m)} ${U(1)} por conglomerado e um ICC de ${fmt(W.icc,2)} (Kish, 1965). `:''}Considerando uma taxa de resposta esperada de ${fmt(W.resp)}%, serão seleccionados ${fmt(R.field)} ${U(1)}${R.clusters?` em ${fmt(R.clusters)} conglomerados`:''}.</p>`;
    if(strat) calc+=`<p>A amostra foi distribuída pelos estratos de forma ${{prop:'proporcional à população',raiz:'proporcional à raiz quadrada da população (Kish, 1965)',igual:'igual',neyman:'óptima de Neyman (Neyman, 1934)'}[ST.alloc]}${ST.target==='estrato'&&prop?', com o objectivo de obter a margem indicada em cada estrato':''}${ST.min>0?`, com um mínimo de ${fmt(ST.min)} ${U(1)} por estrato`:''}. ${prop?`A margem de erro esperada é de ±${fmt(R.overall*100,1)} pontos percentuais para o total e varia entre ±${fmt(Math.min(...R.rows.map(r=>r.marg))*100,1)} e ±${fmt(Math.max(...R.rows.map(r=>r.marg))*100,1)} nos estratos.`:''} A alocação detalhada consta do Anexo.</p>`;
    sec.push(H('4. Tamanho da amostra')+calc);
    sec.push(H('5. Procedimento de selecção')+`<p>${W.design==='aas'?`As unidades da lista serão numeradas de 1 a N e serão sorteados ${fmt(R.field)} números sem reposição, com um gerador pseudo-aleatório e a semente ${W.seed}, o que permite reproduzir o sorteio.`:W.design==='sist'?`Será aplicado um intervalo de amostragem k = N/n, com início aleatório gerado com a semente ${W.seed}.`:W.design==='estrat'?`Em cada estrato será feita uma selecção aleatória independente, com sementes registadas.`:`Na primeira etapa, ${W.design==='multi'?'as áreas serão seleccionadas com probabilidade proporcional ao tamanho (PPS sistemática)':'os conglomerados serão sorteados por aleatória simples'}${strat?', de forma independente em cada estrato':''}, com a semente ${W.seed}. Na segunda etapa, serão seleccionados ${fmt(W.m)} ${U(1)} por amostragem sistemática a partir da listagem actualizada (ICF International, 2012).`} Não haverá substituição de unidades que não respondam.</p>`);
    sec.push(H('6. Ponderação e análise')+`<p>${W.design==='aas'||W.design==='sist'?'Como todas as unidades têm a mesma probabilidade de selecção, a amostra é autoponderada.':`Os pesos de desenho serão calculados como o inverso da probabilidade de selecção${W.design==='multi'?', w = 1/(π₁·π₂)':''}, ajustados pela não resposta. A análise terá em conta os estratos, os conglomerados e os pesos, usando métodos para amostras complexas (Lumley, 2010; Valliant et al., 2018).`} As estimativas serão apresentadas com intervalos de confiança de ${cf}%.</p>`);
  } else {
    const xe=x.en;
    const T=(v,ph)=>v&&v.trim()?esc(v.trim()):`<mark>[${ph}]</mark>`;
    sec.push(H('1. Study population')+`<p>${T(W.title,'study objective')}</p><p>The target population is: ${W.pop.trim()?esc(lc(W.pop)):'<mark>[define the target population]</mark>'} The unit of analysis is the ${U(2)}.</p><p><b>Inclusion criteria:</b> ${T(W.incl,'define')}</p><p><b>Exclusion criteria:</b> ${T(W.excl,'define')}</p>${R.N>0?`<p>The eligible population was estimated at ${fmtEn(R.N)} ${U(3)}${W.nmode==='censo'?`, applying ${fmtEn(W.npct,1)}% to the 2017 Census population (INE)`:''}.</p>`:''}`);
    sec.push(H('2. Sampling frame')+`<p>${W.frame==='lista'?'The sampling frame is the complete list of units in the target population, checked for duplicates and omissions before selection.':W.frame==='areas'?'The first-stage sampling frame is the list of areas (or other clusters) with their size. An updated listing of units will be carried out in the selected areas and used as the second-stage frame (Groves et al., 2009).':'No prior list of units exists. A list of areas (for example, census enumeration areas) will be used as the first-stage frame, followed by a full listing of units in the selected areas (United Nations, 2005).'}</p>`);
    sec.push(H('3. Sampling design')+`<p>A probability sample will be drawn using ${xe}${strat&&W.design==='multi'?', stratified by '+lvl:''}: ${x.defEn}.${strat?` Strata are ${R.rows.length} ${lvl}${ST.level==='dist'&&ST.prov!=='*'?' of '+ST.prov+' province':''}.`:''}</p>`);
    let calc=`<p>Sample size was calculated to estimate ${prop?`an expected proportion of ${fmtEn(W.p)}% with a margin of error of ±${fmtEn(W.d,1)} percentage points`:`a mean with a margin of error of ±${fmtEn(W.dMean,2)} units, assuming a standard deviation of ${fmtEn(W.sigma,2)}`} at ${cf}% confidence (Cochran, 1977), giving a base sample of ${Math.ceil(R.n0)}. ${(prop?W.src:W.src2)?`The expected value is based on: ${esc(prop?W.src:W.src2)}. `:''}${R.N>0?`A finite population correction was applied (N = ${fmtEn(R.N)}), giving ${Math.ceil(R.nF)}. `:''}${isClu()?`To account for clustering, a design effect of ${fmtEn(R.deff,2)} was applied, computed as 1 + (m − 1)·ICC with ${fmtEn(W.m)} ${U(3)} per cluster and an ICC of ${fmtEn(W.icc,2)} (Kish, 1965). `:''}Allowing for an expected response rate of ${fmtEn(W.resp)}%, ${fmtEn(R.field)} ${U(3)} will be selected${R.clusters?` in ${fmtEn(R.clusters)} clusters`:''}.</p>`;
    if(strat) calc+=`<p>The sample was allocated across strata ${{prop:'proportionally to population',raiz:'proportionally to the square root of population (Kish, 1965)',igual:'equally',neyman:'using Neyman optimal allocation (Neyman, 1934)'}[ST.alloc]}${ST.min>0?`, with a minimum of ${fmtEn(ST.min)} ${U(3)} per stratum`:''}. ${prop?`The expected margin of error is ±${fmtEn(R.overall*100,1)} percentage points overall, ranging from ±${fmtEn(Math.min(...R.rows.map(r=>r.marg))*100,1)} to ±${fmtEn(Math.max(...R.rows.map(r=>r.marg))*100,1)} across strata.`:''} The detailed allocation is given in the Annex.</p>`;
    sec.push(H('4. Sample size')+calc);
    sec.push(H('5. Selection procedure')+`<p>${W.design==='aas'?`Units will be numbered 1 to N and ${fmtEn(R.field)} numbers drawn without replacement using a pseudo-random generator with seed ${W.seed}, so the draw can be reproduced.`:W.design==='sist'?`A sampling interval k = N/n will be applied with a random start generated with seed ${W.seed}.`:W.design==='estrat'?'An independent random selection will be made within each stratum, with recorded seeds.':`At the first stage, ${W.design==='multi'?'areas will be selected with probability proportional to size (systematic PPS)':'clusters will be selected by simple random sampling'}${strat?', independently within each stratum':''}, using seed ${W.seed}. At the second stage, ${fmtEn(W.m)} ${U(3)} will be selected by systematic sampling from the updated listing (ICF International, 2012).`} Non-responding units will not be replaced.</p>`);
    sec.push(H('6. Weighting and analysis')+`<p>${W.design==='aas'||W.design==='sist'?'As all units have the same selection probability, the sample is self-weighting.':`Design weights will be computed as the inverse of the selection probability${W.design==='multi'?', w = 1/(π₁·π₂)':''}, adjusted for non-response. Analysis will account for strata, clusters and weights using complex survey methods (Lumley, 2010; Valliant et al., 2018).`} Estimates will be reported with ${cf}% confidence intervals.</p>`);
  }
  return sec.join('');
}
function refsFor(){const k=['groves','cochran'];if(isClu())k.push('kish');if(isStrat())k.push('kish','neyman');if(W.design==='multi')k.push('un2005','dhs2012','valliant');if(W.design!=='aas'&&W.design!=='sist')k.push('lumley');k.push('strobe');return [...new Set(k)];}
function checklist(R){
  const prop=W.out==='prop';
  const items=[
    [W.pop.trim().length>20,'População-alvo definida com pessoa, lugar e tempo'],
    [W.incl.trim().length>10,'Critérios de inclusão escritos'],
    [W.excl.trim().length>10,'Critérios de exclusão escritos'],
    [!!W.frame,'Base de amostragem descrita, com as suas falhas de cobertura'],
    [!!W.design,'Desenho amostral probabilístico escolhido e justificado'],
    [(prop?W.src:W.src2).trim().length>5,'Fonte do valor esperado (p ou σ) indicada'],
    [W.resp<100,'Não resposta prevista no cálculo'],
    [!isClu()||W.icc>0,'Efeito de desenho justificado (ICC e unidades por conglomerado)'],
    [!!W.seed,'Procedimento de selecção reprodutível (semente registada)'],
    [true,'Plano de ponderação e de análise para o desenho escolhido'],
    [false,'Revisão por um estatístico antes da submissão']];
  return items.map(([ok,t])=>`<li class="${ok?'ok':'todo'}"><span aria-hidden="true">${ok?'✓':'!'}</span>${t}${ok?'':' <em>(por fazer)</em>'}</li>`).join('');
}

/* ---------------- summary ---------------- */
function renderSummary(R){
  R=R||calc();
  const row=(k,v,ok)=>`<div class="srow"><dt>${k}</dt><dd class="${ok===false?'miss':''}">${v}</dd></div>`;
  const N=getN();
  $('wz-sum').innerHTML=`<div class="sec-t">O seu estudo</div><dl class="sum">
   ${row('População',W.pop.trim()?esc(W.pop.length>90?W.pop.slice(0,90)+'…':W.pop):'por definir',!!W.pop.trim())}
   ${row('N elegível',N>0?fmt(N):'desconhecido')}
   ${row('Inclusão / exclusão',(W.incl?'✓':'—')+' / '+(W.excl?'✓':'—'),!!(W.incl&&W.excl))}
   ${row('Base de amostragem',{lista:'Lista completa',areas:'Lista de áreas',nenhuma:'Sem lista'}[W.frame])}
   ${row('Indicador',W.out==='prop'?`Proporção, ${fmt(W.p)}% ± ${fmt(W.d,1)} p.p.`:`Média, ± ${fmt(W.dMean,2)}`)}
   ${row('Desenho',D[W.design].n+(W.design==='multi'&&W.stratify?' estratificada':''))}
   ${isStrat()?row('Estratos',`${strataUnits().length} ${ST.level==='prov'?'províncias':'distritos'}`):''}
   </dl>
   <div class="sum-n">${R.error?'<span class="note">'+R.error+'</span>':`<b>${fmt(R.field)}</b><span>${U(1)} a seleccionar${R.clusters?`<br>em ${fmt(R.clusters)} conglomerados`:''}</span>`}</div>
   ${R.error?'':`<div class="note">${fmt(R.days)} dias de campo · ${W.cost>0?fmt(R.cost)+' MT':'custo por definir'}</div>`}`;
}

/* ---------------- navigation ---------------- */
function go(n){
  W.step=Math.max(1,Math.min(6,n));
  document.querySelectorAll('.wz-step').forEach(s=>s.hidden=+s.dataset.step!==W.step);
  document.querySelectorAll('#wz-steps button').forEach(b=>{const k=+b.dataset.go;b.setAttribute('aria-current',k===W.step?'step':'false');b.classList.toggle('done',k<W.step);});
  $('wz-prev').hidden=W.step===1; $('wz-next').textContent=W.step===6?'Voltar ao início':'Seguinte'; $('wz-pos').textContent=`Passo ${W.step} de 6`;
  if(W.step===3) renderDesigns();
  update();
  const top=$('assistente').getBoundingClientRect().top; if(top<-40) $('assistente').scrollIntoView({behavior:'smooth',block:'start'});
}
function update(){
  const R=calc();
  renderParams(R);
  if(W.step===5) renderResult(R);
  if(W.step===6){
    $('w-proto').innerHTML=protocol(R,W.lang);
    document.querySelectorAll('[data-wl]').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.wl===W.lang)));
    $('w-check').innerHTML=checklist(R);
    $('w-refs').innerHTML=refsFor().map(k=>`<li>${REFS[k].f}</li>`).join('');
  }
  const N=getN(); $('w-n-out').innerHTML=W.nmode==='censo'?(isStrat()?`N elegível nos estratos escolhidos: <b>${fmt(N)}</b> (${fmt(W.npct,1)}% da população de 2017 dos estratos).`:`N elegível estimado: <b>${fmt(N)}</b> = ${fmt(areaPop(W.narea))} × ${fmt(W.npct,1)}%.`):W.nmode==='grande'?'Sem correcção para população finita. Fica do lado seguro: a amostra só pode diminuir com N conhecido.':'';
  renderSummary(R);
}

function init(){
  bindText('w-title','title');bindText('w-pop','pop');bindText('w-incl','incl');bindText('w-excl','excl');bindText('w-src','src');bindText('w-src2','src2');
  $('w-unit').value=W.unit;$('w-unit').addEventListener('change',e=>{W.unit=e.target.value;update();});
  $('w-n-area').innerHTML=areaOptions();$('w-n-area').value=W.narea;$('w-n-area').addEventListener('change',e=>{W.narea=e.target.value;update();});
  bindNum('w-n-pct','npct');
  $('w-n-val').value=W.nval;$('w-n-val').addEventListener('input',e=>{W.nval=e.target.value;update();});
  bindChoice('w-nmode','nmode');bindChoice('w-frame','frame');bindChoice('w-obj','obj');bindChoice('w-out','out');bindChoice('w-geo','geo');bindChoice('w-sub','sub');
  bindNum('w-p','p','w-p-r');bindNum('w-sigma','sigma');bindNum('w-resp','resp','w-resp-r');bindNum('w-m','m');bindNum('w-icc','icc','w-icc-r');
  bindNum('w-cost','cost');bindNum('w-perday','perday');bindNum('w-team','team');
  const dEl=$('w-d'),dR=$('w-d-r');
  dEl.addEventListener('input',()=>{const v=parseFloat(dEl.value);if(W.out==='prop'){W.d=v;dR.value=v;}else W.dMean=v;update();});
  dR.addEventListener('input',()=>{W.d=parseFloat(dR.value);dEl.value=dR.value;update();});
  $('w-conf').value=String(W.conf);$('w-conf').addEventListener('change',e=>{W.conf=+e.target.value;update();});
  $('w-stratify').addEventListener('change',e=>{W.stratify=e.target.checked;update();});
  document.querySelectorAll('#wz-steps [data-go]').forEach(b=>b.addEventListener('click',()=>go(+b.dataset.go)));
  $('wz-prev').addEventListener('click',()=>go(W.step-1));
  $('wz-next').addEventListener('click',()=>go(W.step===6?1:W.step+1));
  document.querySelectorAll('[data-wl]').forEach(b=>b.addEventListener('click',()=>{W.lang=b.dataset.wl;update();}));
  $('w-copy').addEventListener('click',e=>copyText(e.currentTarget,()=>$('w-proto').innerText,$('w-proto')));
  $('w-dl').addEventListener('click',()=>download('amostralab-seccao-amostragem.txt',$('w-proto').innerText));
  $('wz-reset').addEventListener('click',()=>{Object.assign(W,JSON.parse(JSON.stringify(EX)));['w-unit'].forEach(id=>$(id).value=W.unit);['w-title','w-pop','w-incl','w-excl','w-src','w-src2'].forEach(id=>$(id).value='');syncChoices();go(1);});
  initStrataUI(); syncChoices(); go(1);
}
A.addStratum=(key,prov)=>{
  if(W.design!=='estrat'&&!(W.design==='multi'&&W.stratify)){W.design='estrat';}
  if(ST.level!=='dist'){ST.level='dist';ST.prov=prov;ST.sel=new Set();}
  else if(ST.prov!=='*'&&ST.prov!==prov){ST.prov='*';}
  ST.sel.add(key); syncStrata(); renderUnits(); go(4);
};
A.addProvince=(prov)=>{if(W.design!=='estrat'&&!(W.design==='multi'&&W.stratify))W.design='estrat';ST.level='dist';ST.prov=prov;selectAllInView();syncStrata();renderUnits();go(4);};
function start(d){DATA=d;init();}
if(A.DATA) start(A.DATA); else A.onData=start;
})();
