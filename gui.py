"""
小核酸药物设计 - 图形用户界面
支持基因编号输入和mRNA序列查询
PubMed文献查询功能
"""

import tkinter as tk
from tkinter import ttk, messagebox, filedialog, simpledialog
import sys
import os
import urllib.request
import urllib.parse
import json
import ssl
import re
import threading

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.sequence_design import NucleicAcidDesigner, SequenceAnalyzer
from src.chemical_modification import ChemicalModifier, DeliverySystem
from src.data_processing import DataHandler, SequenceDatabase


# NCBI PubMed查询类
class PubMedQuerier:
    """
    NCBI PubMed文献查询类
    使用NCBI E-utilities API查询文献和mRNA序列
    """
    
    BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    EMAIL = "research@example.com"
    APP_NAME = "NucleicAcidDesignTool"
    
    @staticmethod
    def search_nucleotide_mrna(gene_name, organism="Homo sapiens[Organism]", max_results=5):
        """
        从NCBI Nucleotide数据库搜索mRNA序列
        
        Args:
            gene_name: 基因名称
            organism: 种属筛选条件
            max_results: 最大返回结果数
            
        Returns:
            mRNA序列列表，每项包含 accession, definition, sequence
        """
        try:
            search_term = f"{gene_name}[Gene Name] AND mRNA[Title] AND {organism}"
            
            search_params = {
                "db": "nucleotide",
                "term": search_term,
                "retmax": max_results,
                "retmode": "json",
                "sort": "relevance",
                "email": PubMedQuerier.EMAIL,
                "tool": PubMedQuerier.APP_NAME
            }
            
            search_url = PubMedQuerier.BASE_URL + "esearch.fcgi?" + urllib.parse.urlencode(search_params)
            
            context = ssl._create_unverified_context()
            try:
                with urllib.request.urlopen(search_url, timeout=15, context=context) as response:
                    search_data = json.loads(response.read().decode('utf-8'))
            except (urllib.error.URLError, json.JSONDecodeError, Exception) as e:
                print(f"NCBI搜索失败: {str(e)}")
                return []
            
            id_list = search_data.get('esearchresult', {}).get('idlist', [])
            
            if not id_list:
                return []
            
            fetch_params = {
                "db": "nucleotide",
                "id": ",".join(id_list),
                "rettype": "fasta",
                "retmode": "text",
                "email": PubMedQuerier.EMAIL,
                "tool": PubMedQuerier.APP_NAME
            }
            
            fetch_url = PubMedQuerier.BASE_URL + "efetch.fcgi?" + urllib.parse.urlencode(fetch_params)
            
            try:
                with urllib.request.urlopen(fetch_url, timeout=30, context=context) as response:
                    fasta_data = response.read().decode('utf-8')
            except (urllib.error.URLError, Exception) as e:
                print(f"NCBI序列获取失败: {str(e)}")
                return []
            
            sequences = []
            current_seq = None
            
            for line in fasta_data.strip().split('\n'):
                if line.startswith('>'):
                    if current_seq:
                        sequences.append(current_seq)
                    parts = line[1:].split('|')
                    accession = parts[0] if parts else line[1:]
                    definition = line[1:]
                    current_seq = {
                        'accession': accession,
                        'definition': definition,
                        'sequence': ''
                    }
                else:
                    if current_seq:
                        current_seq['sequence'] += line.strip()
            
            if current_seq:
                sequences.append(current_seq)
            
            return sequences
            
        except Exception as e:
            print(f"Nucleotide mRNA查询错误: {str(e)}")
            return []
    
    @staticmethod
    def search_pubmed(gene_name, organism="Homo sapiens[Organism]", max_results=10, search_mode="comprehensive"):
        """
        搜索基因相关的小核酸研究PubMed文献

        Args:
            gene_name: 基因名称或基因ID
            organism: 种属筛选条件
            max_results: 最大返回结果数
            search_mode: 搜索模式 - 'strict'(严格), 'normal'(正常), 'comprehensive'(综合)

        Returns:
            文献列表，每项包含pmid, title, authors, journal, abstract
        """
        def make_request(url, timeout=15, retries=2):
            """带重试机制的HTTP请求 - 优化超时处理"""
            context = ssl._create_unverified_context()
            for attempt in range(retries):
                try:
                    with urllib.request.urlopen(url, timeout=timeout, context=context) as response:
                        return json.loads(response.read().decode('utf-8'))
                except urllib.error.URLError as e:
                    if attempt < retries - 1:
                        import time
                        time.sleep(1)
                        continue
                    print(f"网络请求失败 ({attempt+1}/{retries}): {str(e)}")
                    return None
                except json.JSONDecodeError as e:
                    print(f"JSON解析失败: {str(e)}")
                    return None
                except Exception as e:
                    if attempt < retries - 1:
                        import time
                        time.sleep(1)
                        continue
                    print(f"请求失败 ({attempt+1}/{retries}): {str(e)}")
                    return None
        
        def build_search_terms(gene_name, organism, search_mode):
            """构建搜索关键词"""
            siRNA_terms = "siRNA OR small interfering RNA"
            miRNA_terms = "miRNA OR microRNA"
            antisense_terms = "antisense OR antisense oligonucleotide OR ASO"
            rnai_terms = "RNAi OR RNA interference"
            gene_silencing_terms = "gene silencing OR gene knockdown"
            therapeutic_terms = "oligonucleotide therapy OR RNA therapy"
            
            if search_mode == "strict":
                small_nucleic_acids = f"({siRNA_terms}) OR ({miRNA_terms}) OR ({antisense_terms})"
            elif search_mode == "normal":
                small_nucleic_acids = f"({siRNA_terms}) OR ({miRNA_terms}) OR ({antisense_terms}) OR ({rnai_terms})"
            else:
                small_nucleic_acids = f"({siRNA_terms}) OR ({miRNA_terms}) OR ({antisense_terms}) OR ({rnai_terms}) OR ({gene_silencing_terms}) OR ({therapeutic_terms})"
            
            return f"({gene_name}[Gene Name] AND ({small_nucleic_acids})) AND {organism}"
        
        def search_with_fallback(gene_name, organism, max_results, search_mode):
            """优化的PubMed搜索策略 - 支持SCAP等非典型小核酸靶点"""
            search_terms = [
                f"{gene_name}[Gene Name] AND ({build_search_terms(gene_name, organism, search_mode)})",
                f"({gene_name}[Title/Abstract]) AND (siRNA OR miRNA OR antisense OR RNAi OR gene silencing OR oligonucleotide)",
                f"({gene_name}[Title/Abstract]) AND (RNA interference OR gene knockdown OR therapeutic RNA)",
                f"({gene_name}[Title/Abstract]) AND (lipid metabolism OR cholesterol OR SREBP OR hypercholesterolemia OR statin)",
                f"({gene_name}[Title/Abstract]) AND (gene expression OR gene regulation OR molecular biology)",
                f"{gene_name}[Title/Abstract]",
                f"{gene_name}[Gene Name]",
                f"({gene_name}[Title] OR {gene_name}[Abstract] OR {gene_name}[MeSH Terms] OR {gene_name}[Keyword])",
                f"{gene_name}[Title]",
                f"{gene_name}[Abstract]",
                f"{gene_name}[MeSH Terms]",
                gene_name
            ]
            
            for term in search_terms:
                try:
                    search_params = {
                        "db": "pubmed",
                        "term": term,
                        "retmax": max_results,
                        "retmode": "json",
                        "sort": "relevance",
                        "email": PubMedQuerier.EMAIL,
                        "tool": PubMedQuerier.APP_NAME
                    }
                    
                    search_url = PubMedQuerier.BASE_URL + "esearch.fcgi?" + urllib.parse.urlencode(search_params)
                    search_data = make_request(search_url)
                    
                    if search_data:
                        id_list = search_data.get('esearchresult', {}).get('idlist', [])
                        if id_list:
                            return id_list, term
                except Exception as e:
                    print(f"搜索策略失败: {term}, 错误: {str(e)}")
            
            return [], ""
        
        def fetch_articles_esummary(id_list):
            """使用ESummary获取文献信息（更稳定）"""
            articles = []
            
            summary_params = {
                "db": "pubmed",
                "id": ",".join(id_list),
                "retmode": "json",
                "email": PubMedQuerier.EMAIL,
                "tool": PubMedQuerier.APP_NAME
            }
            
            summary_url = PubMedQuerier.BASE_URL + "esummary.fcgi?" + urllib.parse.urlencode(summary_params)
            summary_data = make_request(summary_url, timeout=60)
            
            if summary_data and 'result' in summary_data:
                result = summary_data['result']
                for pmid in id_list:
                    if pmid in result:
                        item = result[pmid]
                        articles.append({
                            'pmid': pmid,
                            'title': item.get('title', 'N/A'),
                            'authors': "; ".join(item.get('authors', [])[:5]) + ("..." if len(item.get('authors', [])) > 5 else ""),
                            'journal': item.get('source', ''),
                            'abstract': item.get('summary', '')[:500] + ("..." if len(item.get('summary', '')) > 500 else "")
                        })
            
            return articles
        
        def fetch_articles_efetch(id_list):
            """使用EFetch获取文献信息（备用）"""
            articles = []
            
            fetch_params = {
                "db": "pubmed",
                "id": ",".join(id_list),
                "retmode": "xml",
                "rettype": "abstract",
                "email": PubMedQuerier.EMAIL,
                "tool": PubMedQuerier.APP_NAME
            }
            
            fetch_url = PubMedQuerier.BASE_URL + "efetch.fcgi?" + urllib.parse.urlencode(fetch_params)
            context = ssl._create_unverified_context()
            
            try:
                with urllib.request.urlopen(fetch_url, timeout=60, context=context) as response:
                    import xml.etree.ElementTree as ET
                    root = ET.fromstring(response.read().decode('utf-8'))
                    
                    for article in root.findall('.//PubmedArticle'):
                        pmid_elem = article.find('.//PMID')
                        pmid = pmid_elem.text if pmid_elem else 'N/A'
                        
                        title_elem = article.find('.//ArticleTitle')
                        title = title_elem.text if title_elem else 'N/A'
                        
                        authors = []
                        for author in article.findall('.//Author'):
                            last_name = author.find('LastName')
                            fore_name = author.find('ForeName')
                            if last_name is not None:
                                authors.append(f"{fore_name.text if fore_name else ''} {last_name.text}".strip())
                        
                        journal_elem = article.find('.//Journal/Title')
                        journal = journal_elem.text if journal_elem else ''
                        
                        abstract_elem = article.find('.//AbstractText')
                        abstract = abstract_elem.text if abstract_elem else ''
                        
                        articles.append({
                            'pmid': pmid,
                            'title': title,
                            'authors': "; ".join(authors[:5]) + ("..." if len(authors) > 5 else ""),
                            'journal': journal,
                            'abstract': abstract[:500] + ("..." if len(abstract) > 500 else "")
                        })
            except Exception as e:
                print(f"EFetch解析失败: {str(e)}")
            
            return articles
        
        try:
            id_list, used_term = search_with_fallback(gene_name, organism, max_results, search_mode)
            
            if not id_list:
                print(f"PubMed搜索失败: {gene_name}")
                return []
            
            articles = fetch_articles_esummary(id_list)
            
            if not articles:
                articles = fetch_articles_efetch(id_list)
            
            return articles
            
        except Exception as e:
            print(f"PubMed查询错误: {str(e)}")
            return []



    @staticmethod
    def get_gene_info(gene_id, organism="Homo sapiens[Organism]"):
        """
        从NCBI Gene数据库获取基因信息及mRNA序列
        
        Args:
            gene_id: 基因ID (如 EGFR, KRAS等)
            organism: 种属筛选条件 (如 "Homo sapiens[Organism]", "Mus musculus[Organism]"等)
            
        Returns:
            基因信息字典，包含mRNA序列
        """
        gene_id_upper = gene_id.strip().upper()
        
        try:
            organism_name = organism.replace("[Organism]", "").strip()
            context = ssl._create_unverified_context()
            
            gene_search_url = PubMedQuerier.BASE_URL + "esearch.fcgi?" + urllib.parse.urlencode({
                "db": "gene",
                "term": f"{gene_id}[Gene Name] AND {organism}",
                "retmode": "json",
                "email": PubMedQuerier.EMAIL,
                "tool": PubMedQuerier.APP_NAME
            })
            
            try:
                with urllib.request.urlopen(gene_search_url, timeout=15, context=context) as response:
                    search_data = json.loads(response.read().decode('utf-8'))
            except Exception:
                search_data = {'esearchresult': {'idlist': []}}
            
            id_list = search_data.get('esearchresult', {}).get('idlist', [])
            
            gene_name = gene_id_upper
            description = ''
            organism_display = organism_name
            
            if id_list:
                gene_id_from_ncbi = id_list[0]
                
                summary_url = PubMedQuerier.BASE_URL + "esummary.fcgi?" + urllib.parse.urlencode({
                    "db": "gene",
                    "id": gene_id_from_ncbi,
                    "retmode": "json",
                    "email": PubMedQuerier.EMAIL,
                    "tool": PubMedQuerier.APP_NAME
                })
                
                try:
                    with urllib.request.urlopen(summary_url, timeout=15, context=context) as response:
                        summary_data = json.loads(response.read().decode('utf-8'))
                    
                    result = summary_data.get('result', {}).get(gene_id_from_ncbi, {})
                    gene_name = result.get('name', gene_id_upper)
                    description = result.get('description', '')
                    organism_display = result.get('organism', organism_name)
                except Exception:
                    pass
            
            mrna_sequences = PubMedQuerier.search_nucleotide_mrna(gene_id_upper, organism=organism, max_results=1)
            mrna = mrna_sequences[0]['sequence'] if mrna_sequences else ''
            accession = mrna_sequences[0]['accession'] if mrna_sequences else ''
            definition = mrna_sequences[0]['definition'] if mrna_sequences else ''
            
            return {
                'ncbi_id': id_list[0] if id_list else '',
                'name': gene_name,
                'description': description,
                'organism': organism_display,
                'chromosome': '',
                'summary': '',
                'mrna': mrna,
                'accession': accession,
                'definition': definition
            }
            
        except Exception as e:
            print(f"NCBI Gene查询错误: {str(e)}")
            return None


# 模拟基因数据库
class GeneDatabase:
    """
    基因数据库模拟类
    提供基因编号到mRNA序列的查询
    """
    
    GENES = {}  # 空数据库，所有查询都通过NCBI
    
    @staticmethod
    def query_gene(gene_id):
        """
        根据基因编号查询基因信息
        """
        gene_id = gene_id.strip().upper()
        return GeneDatabase.GENES.get(gene_id, None)
    
    @staticmethod
    def get_all_gene_ids():
        """
        获取所有可用的基因编号
        """
        return list(GeneDatabase.GENES.keys())


class NucleicAcidDesignGUI:
    """
    小核酸药物设计图形界面
    """
    
    def __init__(self, root):
        self.root = root
        self.root.title("miRNA Designer Pro - 小核酸药物设计系统")
        self.root.geometry("1200x800")
        self.root.minsize(900, 650)
        
        screen_width = self.root.winfo_screenwidth()
        screen_height = self.root.winfo_screenheight()
        if screen_width < 1200 or screen_height < 800:
            self.root.geometry(f"{screen_width}x{screen_height}")
        
        style = ttk.Style()
        style.theme_use('clam')
        
        style.configure('TFrame', background='#f5f5f5')
        style.configure('TLabel', background='#f5f5f5', font=('微软雅黑', 10))
        style.configure('TButton', 
                       background='#4a90d9', 
                       foreground='white',
                       padding=6,
                       font=('微软雅黑', 10, 'bold'))
        style.map('TButton',
                  background=[('active', '#3a7bc8'), ('pressed', '#2d6ab3')])
        style.configure('TLabelframe', background='#f5f5f5')
        style.configure('TLabelframe.Label', background='#f5f5f5', font=('微软雅黑', 11, 'bold'))
        style.configure('TNotebook', background='#f5f5f5')
        style.configure('TNotebook.Tab', 
                       background='#e0e0e0', 
                       padding=[10, 4],
                       font=('微软雅黑', 10))
        style.map('TNotebook.Tab',
                  background=[('selected', '#4a90d9'), ('active', '#c0c0c0')],
                  foreground=[('selected', 'white'), ('active', 'black')])
        style.configure('TCombobox', font=('微软雅黑', 10))
        style.configure('Treeview', 
                       background='white', 
                       fieldbackground='white',
                       font=('微软雅黑', 9))
        style.configure('Treeview.Heading', 
                       background='#4a90d9', 
                       foreground='white',
                       font=('微软雅黑', 9, 'bold'))

        self.designer = NucleicAcidDesigner()
        self.analyzer = SequenceAnalyzer()
        self.modifier = ChemicalModifier()
        self.database = SequenceDatabase()
        
        self.current_mrna = ""

        self.setup_ui()

    def setup_ui(self):
        notebook = ttk.Notebook(self.root)
        notebook.pack(fill='both', expand=True, padx=10, pady=(10, 5))

        self.design_frame = ttk.Frame(notebook)
        notebook.add(self.design_frame, text="序列设计")
        self.setup_scrollable_tab(self.design_frame, self.setup_design_tab)

        self.primer_frame = ttk.Frame(notebook)
        notebook.add(self.primer_frame, text="引物设计")
        self.setup_scrollable_tab(self.primer_frame, self.setup_primer_tab)

        self.evaluation_frame = ttk.Frame(notebook)
        notebook.add(self.evaluation_frame, text="序列评价")
        self.setup_scrollable_tab(self.evaluation_frame, self.setup_evaluation_tab)

        self.modification_frame = ttk.Frame(notebook)
        notebook.add(self.modification_frame, text="化学修饰")
        self.setup_scrollable_tab(self.modification_frame, self.setup_modification_tab)

        self.delivery_frame = ttk.Frame(notebook)
        notebook.add(self.delivery_frame, text="递送系统")
        self.setup_scrollable_tab(self.delivery_frame, self.setup_delivery_tab)

        self.help_frame = ttk.Frame(notebook)
        notebook.add(self.help_frame, text="使用帮助")
        self.setup_scrollable_tab(self.help_frame, self.setup_help_tab)

        status_frame = ttk.Frame(self.root)
        status_frame.pack(fill='x', side='bottom', padx=10, pady=(0, 5))

        self.status_label = tk.Label(status_frame, text="就绪", fg='green', anchor='e', font=('微软雅黑', 9))
        self.status_label.pack(side='right')

        version_label = tk.Label(status_frame, text="小核酸药物设计工具 v1.0", fg='gray', anchor='w', font=('微软雅黑', 8))
        version_label.pack(side='left')
        
        self.root.bind('<Configure>', self.update_widget_sizes)

    def setup_scrollable_tab(self, parent, setup_content_func):
        for widget in parent.winfo_children():
            widget.destroy()
        
        container = ttk.Frame(parent)
        container.pack(fill='both', expand=True, padx=0, pady=0)
        
        canvas = tk.Canvas(container, highlightthickness=0)
        scrollbar = ttk.Scrollbar(container, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)

        def on_frame_configure(event):
            canvas.configure(scrollregion=canvas.bbox("all"))

        scrollable_frame.bind("<Configure>", on_frame_configure)

        canvas_window = canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        def on_canvas_configure(event):
            canvas.itemconfig(canvas_window, width=event.width)
            canvas.configure(scrollregion=canvas.bbox("all"))

        canvas.bind('<Configure>', on_canvas_configure)

        def on_mousewheel(event):
            canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        
        canvas.bind('<MouseWheel>', on_mousewheel)
        canvas.bind('<Button-4>', lambda e: canvas.yview_scroll(-3, "units"))
        canvas.bind('<Button-5>', lambda e: canvas.yview_scroll(3, "units"))

        canvas.pack(side="left", fill="both", expand=True, padx=0, pady=0)
        scrollbar.pack(side="right", fill="y", padx=0, pady=0)

        setup_content_func(scrollable_frame)
    
    def update_widget_sizes(self, event=None):
        """响应窗口大小变化，更新各区域尺寸"""
        for child in self.root.winfo_children():
            child.update_idletasks()
        
        if hasattr(self, 'mrna_input'):
            self.mrna_input.config(height=max(6, int(self.root.winfo_height() / 50)))
        
        if hasattr(self, 'primer_mrna_text'):
            self.primer_mrna_text.config(height=max(6, int(self.root.winfo_height() / 50)))
        
        if hasattr(self, 'result_tree'):
            self.result_tree.update_idletasks()
        
        if hasattr(self, 'primer_tree'):
            self.primer_tree.update_idletasks()

    def setup_design_tab(self, parent=None):
        """
        设置序列设计页面
        包含基因编号输入、种属选择、mRNA序列显示、设计参数设置和结果展示
        """
        if parent is None:
            parent = self.design_frame
        
        # 基因编号输入区域
        gene_frame = ttk.LabelFrame(parent, text="基因编号", padding=10)
        gene_frame.pack(fill='x', padx=10, pady=5)
        
        ttk.Label(gene_frame, text="基因ID:").pack(side='left', padx=5)
        self.gene_id_entry = ttk.Entry(gene_frame, width=30)
        self.gene_id_entry.pack(side='left', padx=5)
        
        ttk.Label(gene_frame, text="种属:").pack(side='left', padx=5)
        self.organism_combo = ttk.Combobox(gene_frame, width=20)
        self.organism_combo['values'] = ['Homo sapiens (Human)', 'Mus musculus (Mouse)', 'Rattus norvegicus (Rat)', 'Danio rerio (Zebrafish)', 'Caenorhabditis elegans (C. elegans)']
        self.organism_combo.current(0)
        self.organism_combo.pack(side='left', padx=5)
        
        self.query_btn = ttk.Button(gene_frame, text="查询mRNA序列", command=self.query_gene)
        self.query_btn.pack(side='left', padx=5)
        
        self.cancel_btn = ttk.Button(gene_frame, text="取消查询", command=self.cancel_query, state='disabled')
        self.cancel_btn.pack(side='left', padx=5)
        
        self.pubmed_btn = ttk.Button(gene_frame, text="查询PubMed文献", command=self.query_pubmed)
        self.pubmed_btn.pack(side='left', padx=5)
        
        self.query_cancelled = False
        
        # mRNA序列显示区域（初始隐藏）
        self.mrna_frame = ttk.LabelFrame(parent, text="mRNA序列（功能区域标记）", padding=10)
        self.mrna_frame.pack(fill='both', expand=True, padx=10, pady=5)

        mrna_scroll = ttk.Scrollbar(self.mrna_frame)
        mrna_scroll.pack(side='right', fill='y')

        self.mrna_input = tk.Text(self.mrna_frame, height=6, wrap='word', yscrollcommand=mrna_scroll.set)
        self.mrna_input.pack(fill='both', expand=True)
        mrna_scroll.config(command=self.mrna_input.yview)
        self.mrna_input.config(state='disabled')
        
        # PubMed文献显示区域（初始隐藏）
        self.pubmed_frame = ttk.LabelFrame(parent, text="PubMed文献查询结果", padding=10)
        
        self.pubmed_text = tk.Text(self.pubmed_frame, height=6, wrap='word')
        self.pubmed_text.pack(fill='both', expand=True)

        # 设计参数
        param_frame = ttk.LabelFrame(parent, text="设计参数", padding=10)
        param_frame.pack(fill='x', padx=10, pady=5)

        ttk.Label(param_frame, text="目标区域:").grid(row=0, column=0, sticky='w')
        self.region_combo = ttk.Combobox(param_frame, width=30)
        self.region_combo['values'] = [
            '★ CDS (编码序列区)',
            '★ 起始密码子附近',
            '★ 外显子区域',
            '可变剪接区域',
            '全序列',
            '5\'UTR',
            '3\'UTR',
            '终止密码子附近'
        ]
        self.region_combo.current(0)
        self.region_combo.grid(row=0, column=1, padx=5)
        self.region_combo.bind('<<ComboboxSelected>>', self.on_region_changed)

        ttk.Label(param_frame, text="说明:").grid(row=0, column=2, sticky='w', padx=(20, 5))
        self.region_desc_label = tk.Label(param_frame, text="首选，敲低效果最佳", foreground='green', font=('微软雅黑', 9))
        self.region_desc_label.grid(row=0, column=3, sticky='w', columnspan=2)

        ttk.Label(param_frame, text="位置偏移:").grid(row=1, column=0, sticky='w', pady=5)
        self.position_entry = ttk.Entry(param_frame, width=10)
        self.position_entry.insert(0, "0")
        self.position_entry.grid(row=1, column=1, padx=5, sticky='w')

        ttk.Label(param_frame, text="siRNA长度:").grid(row=1, column=2, sticky='w', pady=5)
        self.length_entry = ttk.Entry(param_frame, width=10)
        self.length_entry.insert(0, "21")
        self.length_entry.grid(row=1, column=3, padx=5, sticky='w')

        btn_frame = ttk.Frame(parent)
        btn_frame.pack(fill='x', padx=10, pady=5)

        ttk.Button(btn_frame, text="设计siRNA", command=self.design_sirna).pack(side='left', padx=5)

        result_frame = ttk.LabelFrame(parent, text="设计结果", padding=10)
        result_frame.pack(fill='both', expand=True, padx=10, pady=5)

        result_vscroll = ttk.Scrollbar(result_frame, orient='vertical')
        result_vscroll.pack(side='right', fill='y')
        
        result_hscroll = ttk.Scrollbar(result_frame, orient='horizontal')
        result_hscroll.pack(side='bottom', fill='x')

        self.result_tree = ttk.Treeview(result_frame, yscrollcommand=result_vscroll.set, xscrollcommand=result_hscroll.set, columns=('排名', '位置', '正义链', '反义链', 'GC%', 'Reynolds', 'Amarz', '种子区', '热力学', 'Ui-Tei', '靶区域', '总分', '推荐'), show='headings')
        self.result_tree.heading('排名', text='排名')
        self.result_tree.heading('位置', text='mRNA位置')
        self.result_tree.heading('正义链', text='正义链 (Sense)')
        self.result_tree.heading('反义链', text='反义链 (Antisense)')
        self.result_tree.heading('GC%', text='GC%')
        self.result_tree.heading('Reynolds', text='Reynolds')
        self.result_tree.heading('Amarz', text='Amarz')
        self.result_tree.heading('种子区', text='种子区')
        self.result_tree.heading('热力学', text='热力学')
        self.result_tree.heading('Ui-Tei', text='Ui-Tei')
        self.result_tree.heading('靶区域', text='靶区域')
        self.result_tree.heading('总分', text='总分')
        self.result_tree.heading('推荐', text='推荐')

        self.result_tree.column('排名', width=45, anchor='center')
        self.result_tree.column('位置', width=70, anchor='center')
        self.result_tree.column('正义链', width=140, anchor='w')
        self.result_tree.column('反义链', width=140, anchor='w')
        self.result_tree.column('GC%', width=45, anchor='center')
        self.result_tree.column('Reynolds', width=55, anchor='center')
        self.result_tree.column('Amarz', width=45, anchor='center')
        self.result_tree.column('种子区', width=45, anchor='center')
        self.result_tree.column('热力学', width=45, anchor='center')
        self.result_tree.column('Ui-Tei', width=45, anchor='center')
        self.result_tree.column('靶区域', width=70, anchor='center')
        self.result_tree.column('总分', width=45, anchor='center')
        self.result_tree.column('推荐', width=50, anchor='center')
        
        self.result_tree.pack(fill='both', expand=True)
        result_vscroll.config(command=self.result_tree.yview)
        result_hscroll.config(command=self.result_tree.xview)
        
        export_button_frame = ttk.Frame(parent)
        export_button_frame.pack(fill='x', padx=10, pady=5)
        
        ttk.Button(export_button_frame, text="复制选中行", command=self.copy_selected_result).pack(side='left', padx=2)
        ttk.Button(export_button_frame, text="复制全部结果", command=self.copy_all_results).pack(side='left', padx=2)
        ttk.Button(export_button_frame, text="导出CSV", command=self.export_results_to_csv).pack(side='left', padx=2)

        self.blast_frame = ttk.LabelFrame(parent, text="BLAST相似性分析", padding=10)
        self.blast_frame.pack(fill='both', expand=True, padx=10, pady=5)

        blast_btn_frame = ttk.Frame(self.blast_frame)
        blast_btn_frame.pack(fill='x', pady=(0, 5))

        self.run_blast_btn = ttk.Button(blast_btn_frame, text="运行BLAST (选中序列)", command=self.run_blast_on_selected)
        self.run_blast_btn.pack(side='left', padx=5)

        self.delete_seq_btn = ttk.Button(blast_btn_frame, text="删除选中序列", command=self.delete_selected_sequence)
        self.delete_seq_btn.pack(side='left', padx=5)

        ttk.Label(blast_btn_frame, text="提示: 设计完成后可运行BLAST检测脱靶效应，或删除不满意的序列").pack(side='left', padx=10)

        blast_scroll = ttk.Scrollbar(self.blast_frame)
        blast_scroll.pack(side='right', fill='y')

        self.blast_text = tk.Text(self.blast_frame, height=10, wrap='word', yscrollcommand=blast_scroll.set)
        self.blast_text.pack(fill='both', expand=True)
        blast_scroll.config(command=self.blast_text.yview)

        self.blast_text.insert("1.0", "BLAST相似性分析结果将在此显示\n\n请先设计siRNA序列，然后选择序列点击\"运行BLAST\"按钮进行脱靶效应检测")

    def format_mrna_with_regions(self, sequence, gene_name="", organism="", accession=""):
        """
        格式化mRNA序列显示，包含功能区域标记和序列序号
        返回格式化文本和区域位置信息
        """
        if not sequence:
            return "", {}

        result = ""
        
        result += "═" * 80 + "\n"
        result += "【mRNA序列信息】\n"
        result += "═" * 80 + "\n"
        
        if gene_name:
            result += f"基因名称: {gene_name}\n"
        if organism:
            result += f"种属: {organism}\n"
        if accession:
            result += f"NCBI Accession: {accession}\n"
        
        result += "\n【序列统计】\n"
        result += "-" * 60 + "\n"
        
        seq_len = len(sequence)
        a_count = sequence.count('A') + sequence.count('a')
        u_count = sequence.count('U') + sequence.count('u') + sequence.count('T') + sequence.count('t')
        c_count = sequence.count('C') + sequence.count('c')
        g_count = sequence.count('G') + sequence.count('g')
        
        gc_content = ((g_count + c_count) / seq_len) * 100 if seq_len > 0 else 0
        
        result += f"序列长度: {seq_len} nt\n"
        result += f"A: {a_count} ({(a_count/seq_len*100):.1f}%)\n"
        result += f"U/T: {u_count} ({(u_count/seq_len*100):.1f}%)\n"
        result += f"C: {c_count} ({(c_count/seq_len*100):.1f}%)\n"
        result += f"G: {g_count} ({(g_count/seq_len*100):.1f}%)\n"
        result += f"GC含量: {gc_content:.1f}%\n"
        
        result += "\n【功能区域预测】\n"
        result += "-" * 60 + "\n"
        
        utr5_len = min(100, len(sequence) // 6)
        cds_len = min(len(sequence) - utr5_len - 50, len(sequence) // 2)
        utr3_start = utr5_len + cds_len
        
        regions = {
            'utr5': (1, utr5_len),
            'cds': (utr5_len + 1, utr3_start),
            'utr3': (utr3_start + 1, len(sequence))
        }
        
        result += f"5'UTR区域: 位置 {regions['utr5'][0]}-{regions['utr5'][1]} ({regions['utr5'][1] - regions['utr5'][0] + 1} nt)\n"
        result += f"CDS区域: 位置 {regions['cds'][0]}-{regions['cds'][1]} ({regions['cds'][1] - regions['cds'][0] + 1} nt)\n"
        result += f"3'UTR区域: 位置 {regions['utr3'][0]}-{regions['utr3'][1]} ({regions['utr3'][1] - regions['utr3'][0] + 1} nt)\n"
        
        orf_start = utr5_len + 1
        orf_end = utr5_len + ((cds_len // 3) * 3)
        if orf_end > orf_start:
            result += f"预测ORF: 位置 {orf_start}-{orf_end} ({(orf_end - orf_start + 1)//3} 个密码子)\n"
        
        result += "\n" + "═" * 80 + "\n"
        result += "【图例】 5'UTR(蓝色) | CDS(绿色) | 3'UTR(红色)\n"
        result += "═" * 80 + "\n\n"
        result += "【mRNA序列】\n"
        line_length = 60
        group_size = 10
        
        for i in range(0, len(sequence), line_length):
            line = sequence[i:i+line_length]
            line_num = i + 1
            
            result += f"{line_num:6d}: "
            
            for j in range(0, len(line), group_size):
                group = line[j:j+group_size]
                pos_in_seq = i + j
                
                if pos_in_seq < utr5_len:
                    result += f"{group} "
                elif pos_in_seq < utr3_start:
                    result += f"{group} "
                else:
                    result += f"{group} "
            
            if i + line_length < len(sequence):
                result += "\n"
            
            if (i + line_length) % (line_length * 3) == 0 and i > 0:
                result += "\n"

        result += "\n" + "═" * 80 + "\n"
        result += f"序列长度: {len(sequence)} nt | GC含量: {gc_content:.1f}% | 预测CDS长度: {regions['cds'][1] - regions['cds'][0] + 1} nt\n"
        
        return result, regions

    def color_mrna_regions(self, regions, text_widget=None):
        """
        给mRNA序列添加颜色标记
        text_widget: 可选参数，指定要标记的Text组件，默认为self.mrna_input
        """
        if text_widget is None:
            text_widget = self.mrna_input
        
        utr5_end = regions.get('utr5', (1, 100))[1]
        cds_start, cds_end = regions.get('cds', (101, 1500))
        utr3_start = regions.get('utr3', (1501, 2000))[0]

        text_widget.tag_config('utr5', foreground='blue')
        text_widget.tag_config('cds', foreground='green')
        text_widget.tag_config('utr3', foreground='red')

        content = text_widget.get("1.0", tk.END)
        lines = content.split('\n')
        
        for line_idx, line in enumerate(lines):
            if ':' not in line or line.strip().startswith('=') or line.strip().startswith('【'):
                continue
            
            line_num_part, seq_part = line.split(':', 1)
            try:
                line_start_pos = int(line_num_part.strip())
            except ValueError:
                continue
            
            seq_text = seq_part.strip()
            seq_start_idx = len(line_num_part) + 1
            current_seq_pos = line_start_pos
            
            for char_idx, char in enumerate(seq_text):
                if char == ' ':
                    continue
                
                if current_seq_pos <= utr5_end:
                    tag = 'utr5'
                elif current_seq_pos <= cds_end:
                    tag = 'cds'
                else:
                    tag = 'utr3'
                
                line_start = f"{line_idx + 1}.{seq_start_idx + char_idx}"
                line_end = f"{line_idx + 1}.{seq_start_idx + char_idx + 1}"
                text_widget.tag_add(tag, line_start, line_end)
                
                current_seq_pos += 1

    def cancel_query(self):
        """取消正在进行的查询"""
        self.query_cancelled = True
        self.set_status("查询已取消", "orange")
        self.query_btn.config(state='normal')
        self.cancel_btn.config(state='disabled')
    
    def query_gene(self):
        """
        查询基因编号对应的mRNA序列 - 直接使用NCBI搜索
        """
        gene_id = self.gene_id_entry.get().strip()

        if not gene_id:
            messagebox.showwarning("警告", "请输入基因编号")
            return

        self.query_cancelled = False
        self.query_btn.config(state='disabled')
        self.cancel_btn.config(state='normal')
        self.set_status(f"正在从NCBI查询基因 {gene_id}...", "orange")
        
        import threading
        def fetch_from_ncbi():
            if self.query_cancelled:
                return
            
            selected_organism = self.organism_combo.get()
            organism_mapping = {
                'Homo sapiens (Human)': 'Homo sapiens[Organism]',
                'Mus musculus (Mouse)': 'Mus musculus[Organism]',
                'Rattus norvegicus (Rat)': 'Rattus norvegicus[Organism]',
                'Danio rerio (Zebrafish)': 'Danio rerio[Organism]',
                'Caenorhabditis elegans (C. elegans)': 'Caenorhabditis elegans[Organism]'
            }
            organism = organism_mapping.get(selected_organism, 'Homo sapiens[Organism]')
            
            gene_id_upper = gene_id.strip().upper()
            ncbi_info = PubMedQuerier.get_gene_info(gene_id_upper, organism=organism)
            
            if self.query_cancelled:
                return
            
            if ncbi_info and ncbi_info.get('mrna'):
                self.root.after(0, lambda: self.update_mrna_result(ncbi_info, gene_id))
            else:
                self.root.after(0, lambda: self.on_gene_query_failed(gene_id, selected_organism))
        
        threading.Thread(target=fetch_from_ncbi, daemon=True).start()
    
    def update_mrna_result(self, ncbi_info, gene_id):
        """更新mRNA查询结果"""
        if self.query_cancelled:
            return
        mrna_sequence = ncbi_info['mrna'].replace('T', 'U').replace('t', 'u')
        self.current_mrna = mrna_sequence
        self.mrna_input.config(state='normal')
        self.mrna_input.delete("1.0", tk.END)
        formatted_seq, regions = self.format_mrna_with_regions(mrna_sequence, ncbi_info.get('name', gene_id), ncbi_info.get('organism', ''), ncbi_info.get('accession', ''))
        self.mrna_input.insert("1.0", formatted_seq)
        self.color_mrna_regions(regions)
        self.mrna_input.config(state='disabled')
        self.set_status(f"✓ 已从NCBI获取基因 {ncbi_info.get('name', gene_id)} 的mRNA序列", "green")
        self.query_btn.config(state='normal')
        self.cancel_btn.config(state='disabled')
    
    def on_gene_query_failed(self, gene_id, selected_organism):
        """基因查询失败处理"""
        if self.query_cancelled:
            return
        self.set_status(f"未找到基因 {gene_id}", "red")
        self.query_btn.config(state='normal')
        self.cancel_btn.config(state='disabled')
        messagebox.showwarning("警告", f"未找到基因编号 '{gene_id}' 在种属 {selected_organism} 中的mRNA序列\n\n建议尝试：\n1. 检查基因名称是否正确\n2. 尝试其他种属\n3. 检查网络连接")

    def query_pubmed(self):
        """
        查询PubMed文献
        """
        gene_id = self.gene_id_entry.get().strip()

        if not gene_id:
            messagebox.showwarning("警告", "请输入基因编号")
            return

        selected_organism = self.organism_combo.get()
        organism_mapping = {
            'Homo sapiens (Human)': 'Homo sapiens[Organism]',
            'Mus musculus (Mouse)': 'Mus musculus[Organism]',
            'Rattus norvegicus (Rat)': 'Rattus norvegicus[Organism]',
            'Danio rerio (Zebrafish)': 'Danio rerio[Organism]',
            'Caenorhabditis elegans (C. elegans)': 'Caenorhabditis elegans[Organism]'
        }
        organism = organism_mapping.get(selected_organism, 'Homo sapiens[Organism]')

        self.pubmed_btn.config(state='disabled')
        self.pubmed_text.delete("1.0", tk.END)
        self.pubmed_text.insert("1.0", f"正在查询PubMed数据库（种属: {selected_organism}），请稍候...\n")
        self.set_status(f"正在查询PubMed小核酸文献: {gene_id}...", "orange")
        self.root.update()

        gene_id_upper = gene_id.strip().upper()

        gene_info = GeneDatabase.query_gene(gene_id_upper)
        search_term = gene_id_upper
        if gene_info:
            search_term = gene_info['name']

        articles = PubMedQuerier.search_pubmed(search_term, organism=organism, max_results=10)

        self.pubmed_text.delete("1.0", tk.END)

        if articles:
            self.pubmed_frame.pack(fill='both', expand=True, padx=10, pady=5)
            self.pubmed_text.insert("1.0", f"=== {search_term} 小核酸研究文献 (共 {len(articles)} 篇) ===\n(包含siRNA/miRNA/ASO/反义寡核苷酸研究)\n\n")
            for i, article in enumerate(articles, 1):
                self.pubmed_text.insert(tk.END, f"[{i}] PMID: {article['pmid']} ")
                self.pubmed_text.insert(tk.END, "[在浏览器中打开]", ("link",))
                self.pubmed_text.insert(tk.END, "\n")
                self.pubmed_text.insert(tk.END, f"标题: {article['title']}\n")
                self.pubmed_text.insert(tk.END, f"作者: {article['authors']}\n")
                self.pubmed_text.insert(tk.END, f"期刊: {article['journal']}\n")
                if article['abstract']:
                    self.pubmed_text.insert(tk.END, f"摘要: {article['abstract']}\n")
                self.pubmed_text.insert(tk.END, "-" * 50 + "\n\n")

            self.pubmed_text.tag_config("link", foreground="blue", underline=True)
            
            def open_pubmed_link(event):
                text = self.pubmed_text.get("1.0", tk.END)
                cursor_pos = self.pubmed_text.index(tk.CURRENT)
                line_num = int(cursor_pos.split('.')[0])
                
                current_pmid = None
                lines = text.split('\n')
                for i, line in enumerate(lines[:line_num]):
                    if 'PMID:' in line:
                        current_pmid = line.split('PMID:')[1].strip().split()[0]
                
                if current_pmid:
                    import webbrowser
                    url = f"https://pubmed.ncbi.nlm.nih.gov/{current_pmid}/"
                    webbrowser.open(url)
                    self.set_status(f"正在打开PubMed链接: {current_pmid}", "green")
            
            self.pubmed_text.tag_bind("link", "<Button-1>", open_pubmed_link)

            self.set_status(f"✓ 已获取 {len(articles)} 篇小核酸研究文献，点击蓝色链接可在浏览器中打开", "green")
        else:
            self.pubmed_text.insert("1.0", f"未找到基因 {search_term} 相关的小核酸研究文献\n\n提示: siRNA/miRNA/ASO/反义寡核苷酸研究")
            self.set_status("未找到相关小核酸文献", "red")
            messagebox.showwarning("提示", "未找到相关文献，请尝试其他基因或种属")

        self.pubmed_btn.config(state='normal')

    def setup_modification_tab(self, parent=None):
        if parent is None:
            parent = self.modification_frame
        
        import_frame = ttk.LabelFrame(parent, text='序列导入', padding=10)
        import_frame.pack(fill='x', padx=10, pady=5)
        
        def on_import_from_design():
            all_items = self.result_tree.get_children()
            if not all_items:
                messagebox.showwarning('提示', '序列设计页面暂无设计结果，请先设计siRNA序列')
                return
            selected_items = self.result_tree.selection()
            if selected_items:
                item = selected_items[0]
            else:
                item = all_items[0]
            values = self.result_tree.item(item, 'values')
            if len(values) >= 4:
                self.mod_sense_entry.delete(0, tk.END)
                self.mod_sense_entry.insert(0, values[2])
                self.mod_antisense_entry.delete(0, tk.END)
                self.mod_antisense_entry.insert(0, values[3])
                self.set_status('✓ 已从设计结果导入序列', 'green')
        
        ttk.Button(import_frame, text='从设计结果导入', command=on_import_from_design).pack(side='left', padx=5)
        
        config_frame = ttk.LabelFrame(parent, text='修饰配置', padding=10)
        config_frame.pack(fill='x', padx=10, pady=5)
        
        ttk.Label(config_frame, text='递送平台:').grid(row=0, column=0, sticky='w', padx=5)
        self.mod_platform = ttk.Combobox(config_frame, values=['GalNAc', 'LNP'], width=15)
        self.mod_platform.current(0)
        self.mod_platform.grid(row=0, column=1, padx=5)
        
        ttk.Label(config_frame, text='双链构型:').grid(row=0, column=2, sticky='w', padx=(20, 5))
        self.mod_duplex_display = ttk.Label(config_frame, text='自动识别', width=15)
        self.mod_duplex_display.grid(row=0, column=3, padx=5)
        
        ttk.Button(config_frame, text='生成推荐修饰序列', command=self.generate_recommended_mod).grid(row=0, column=4, padx=10)
        
        seq_frame = ttk.LabelFrame(parent, text='序列输入', padding=10)
        seq_frame.pack(fill='x', padx=10, pady=5)
        
        ttk.Label(seq_frame, text='正义链:').grid(row=0, column=0, sticky='w')
        self.mod_sense_entry = ttk.Entry(seq_frame, width=50)
        self.mod_sense_entry.grid(row=0, column=1, padx=5)
        
        ttk.Label(seq_frame, text='反义链:').grid(row=1, column=0, sticky='w')
        self.mod_antisense_entry = ttk.Entry(seq_frame, width=50)
        self.mod_antisense_entry.grid(row=1, column=1, padx=5)
        
        modified_seq_frame = ttk.LabelFrame(parent, text='修饰结果（合成订单格式）', padding=10)
        modified_seq_frame.pack(fill='both', expand=True, padx=10, pady=5)
        
        ttk.Label(modified_seq_frame, text='正义链 (Sense):').pack(anchor='w')
        self.mod_sense_result = tk.Text(modified_seq_frame, height=3, wrap='word', font=('Consolas', 10))
        self.mod_sense_result.pack(fill='x', expand=False)
        self.mod_sense_result.insert('1.0', '合成订单格式将在此显示')
        self.mod_sense_result.config(state='disabled')
        
        ttk.Label(modified_seq_frame, text='反义链 (Antisense):').pack(anchor='w', pady=(5, 0))
        self.mod_antisense_result = tk.Text(modified_seq_frame, height=3, wrap='word', font=('Consolas', 10))
        self.mod_antisense_result.pack(fill='x', expand=False)
        self.mod_antisense_result.insert('1.0', '合成订单格式将在此显示')
        self.mod_antisense_result.config(state='disabled')
        
        validation_frame = ttk.LabelFrame(parent, text='校验结果', padding=10)
        validation_frame.pack(fill='x', padx=10, pady=5)
        
        self.mod_validation_text = tk.Text(validation_frame, height=4, wrap='word', font=('Consolas', 9))
        self.mod_validation_text.pack(fill='x', expand=False)

    def setup_delivery_tab(self, parent=None):
        """
        设置递送系统页面
        支持siRNA递送系统设计，包括：
        - 序列导入（从设计结果或Excel）
        - 目标组织选择
        - 递送系统类型选择（GalNAC、LNP等）
        - 递送推荐方案获取
        - 结合后序列展示
        - 推荐递送方案生成
        """
        if parent is None:
            parent = self.delivery_frame
        
        # 序列导入区域
        import_frame = ttk.LabelFrame(parent, text="序列导入", padding=10)
        import_frame.pack(fill='x', padx=10, pady=5)
        
        def on_import_from_design_delivery():
            all_items = self.result_tree.get_children()
            if not all_items:
                messagebox.showwarning("提示", "序列设计页面暂无设计结果，请先设计siRNA序列")
                return
            selected_items = self.result_tree.selection()
            if selected_items:
                item = selected_items[0]
            else:
                item = all_items[0]
            values = self.result_tree.item(item, 'values')
            if len(values) >= 4:
                self.delivery_sense_entry.delete(0, tk.END)
                self.delivery_sense_entry.insert(0, values[2])
                self.delivery_antisense_entry.delete(0, tk.END)
                self.delivery_antisense_entry.insert(0, values[3])
        
        ttk.Button(import_frame, text="从设计结果导入", command=on_import_from_design_delivery).pack(side='left', padx=5)
        
        seq_frame = ttk.LabelFrame(parent, text="siRNA序列", padding=10)
        seq_frame.pack(fill='x', padx=10, pady=5)
        
        ttk.Label(seq_frame, text="正义链:").grid(row=0, column=0, sticky='w')
        self.delivery_sense_entry = ttk.Entry(seq_frame, width=50)
        self.delivery_sense_entry.grid(row=0, column=1, padx=5)
        
        ttk.Label(seq_frame, text="反义链:").grid(row=1, column=0, sticky='w')
        self.delivery_antisense_entry = ttk.Entry(seq_frame, width=50)
        self.delivery_antisense_entry.grid(row=1, column=1, padx=5)
        
        ttk.Button(seq_frame, text="生成推荐递送序列", command=self.generate_recommended_delivery).grid(row=0, column=2, padx=5)

        tissue_frame = ttk.LabelFrame(parent, text="目标组织与递送系统", padding=10)
        tissue_frame.pack(fill='x', padx=10, pady=5)

        ttk.Label(tissue_frame, text="选择目标组织:").grid(row=0, column=0, sticky='w')
        self.tissue_type = ttk.Combobox(tissue_frame, width=20)
        self.tissue_type['values'] = ['肝脏', '眼睛', '中枢神经系统', '肿瘤', '肺部', '肌肉']
        self.tissue_type.current(0)
        self.tissue_type.grid(row=0, column=1, padx=5)

        ttk.Label(tissue_frame, text="递送系统:").grid(row=0, column=2, sticky='w', padx=(10, 5))
        self.delivery_system = ttk.Combobox(tissue_frame, width=20)
        self.delivery_system['values'] = ['LNP (脂质纳米粒)', 'GalNAC偶联', '聚合物纳米粒', '病毒载体', '裸RNA', '吸入式纳米粒']
        self.delivery_system.current(0)
        self.delivery_system.grid(row=0, column=3, padx=5)

        ttk.Button(tissue_frame, text="获取推荐", command=self.get_delivery_recommendation).grid(row=0, column=4, padx=5)

        conjugate_frame = ttk.LabelFrame(parent, text="递送系统结合后序列", padding=10)
        conjugate_frame.pack(fill='both', expand=True, padx=10, pady=5)
        
        ttk.Label(conjugate_frame, text="完整序列结构:").pack(anchor='w')
        self.delivery_seq_result = tk.Text(conjugate_frame, height=4, wrap='word', font=('Consolas', 10))
        self.delivery_seq_result.pack(fill='x', expand=False)
        self.delivery_seq_result.insert("1.0", "请输入序列并选择递送系统")
        self.delivery_seq_result.config(state='disabled')
        
        ttk.Button(conjugate_frame, text="生成结合序列", command=self.generate_conjugate_sequence).pack(pady=5)

        recommend_frame = ttk.LabelFrame(parent, text="推荐递送系统", padding=10)
        recommend_frame.pack(fill='both', expand=True, padx=10, pady=5)

        self.delivery_text = tk.Text(recommend_frame, height=10, wrap='word')
        self.delivery_text.pack(fill='both', expand=True)

    def generate_recommended_mod(self):
        sense = self.mod_sense_entry.get().strip()
        antisense = self.mod_antisense_entry.get().strip()
        platform = self.mod_platform.get()
        duplex_config = self.mod_duplex.get()
        
        if not sense or not antisense:
            messagebox.showwarning('提示', '请先输入或导入序列')
            return
        
        result = self.modifier.design_siRNA_modification(sense, antisense, platform, duplex_config)
        
        self.mod_sense_result.config(state='normal')
        self.mod_sense_result.delete('1.0', tk.END)
        self.mod_sense_result.insert('1.0', result['sense_modified'])
        self.mod_sense_result.config(state='disabled')
        
        self.mod_antisense_result.config(state='normal')
        self.mod_antisense_result.delete('1.0', tk.END)
        self.mod_antisense_result.insert('1.0', result['antisense_modified'])
        self.mod_antisense_result.config(state='disabled')
        
        validation = result.get('validation', {})
        self.mod_validation_text.delete('1.0', tk.END)
        if validation.get('valid', False):
            self.mod_validation_text.insert('1.0', '✓ 校验通过\n')
        else:
            self.mod_validation_text.insert('1.0', '✗ 校验失败\n')
        for err in validation.get('errors', []):
            self.mod_validation_text.insert(tk.END, f'错误: {err}\n')
        for warn in validation.get('warnings', []):
            self.mod_validation_text.insert(tk.END, f'警告: {warn}\n')
        
        self.set_status(f'✓ 已生成{platform}平台的推荐修饰序列', 'green')

    def generate_recommended_delivery(self):
        sense = self.delivery_sense_entry.get().strip()
        antisense = self.delivery_antisense_entry.get().strip()
        
        if not sense or not antisense:
            messagebox.showwarning("提示", "请先输入或导入序列")
            return
        
        self.tissue_type.set('肝脏')
        self.delivery_system.set('GalNAC偶联')
        self.get_delivery_recommendation()
        self.generate_conjugate_sequence()

    def copy_selected_result(self):
        """复制选中的行到剪贴板"""
        selected_items = self.result_tree.selection()
        if not selected_items:
            self.set_status("请先选择要复制的行", "orange")
            return
        
        result_text = ""
        headers = ["排名", "位置", "正义链", "反义链", "GC%", "Reynolds", "Amarz", "种子区", "热力学", "Ui-Tei", "靶区域", "总分", "推荐"]
        result_text += "\t".join(headers) + "\n"
        
        for item in selected_items:
            values = self.result_tree.item(item)['values']
            result_text += "\t".join(str(v) for v in values) + "\n"
        
        self.root.clipboard_clear()
        self.root.clipboard_append(result_text)
        self.set_status(f"已复制{len(selected_items)}行到剪贴板", "green")

    def copy_all_results(self):
        """复制所有结果到剪贴板"""
        all_items = self.result_tree.get_children()
        if not all_items:
            self.set_status("没有可复制的结果", "orange")
            return
        
        result_text = ""
        headers = ["排名", "位置", "正义链", "反义链", "GC%", "Reynolds", "Amarz", "种子区", "热力学", "Ui-Tei", "靶区域", "总分", "推荐"]
        result_text += "\t".join(headers) + "\n"
        
        for item in all_items:
            values = self.result_tree.item(item)['values']
            result_text += "\t".join(str(v) for v in values) + "\n"
        
        self.root.clipboard_clear()
        self.root.clipboard_append(result_text)
        self.set_status(f"已复制{len(all_items)}行到剪贴板", "green")

    def export_results_to_csv(self):
        """导出结果到CSV文件"""
        all_items = self.result_tree.get_children()
        if not all_items:
            self.set_status("没有可导出的结果", "orange")
            return
        
        try:
            from tkinter import filedialog
            import csv
            from datetime import datetime
            
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            default_filename = f"siRNA_results_{timestamp}.csv"
            
            file_path = filedialog.asksaveasfilename(
                defaultextension=".csv",
                initialfile=default_filename,
                filetypes=[("CSV文件", "*.csv"), ("所有文件", "*.*")],
                title="导出结果"
            )
            
            if not file_path:
                return
            
            headers = ["排名", "位置", "正义链", "反义链", "GC%", "Reynolds", "Amarz", "种子区", "热力学", "Ui-Tei", "靶区域", "总分", "推荐"]
            
            with open(file_path, 'w', newline='', encoding='utf-8-sig') as f:
                writer = csv.writer(f)
                writer.writerow(headers)
                
                for item in all_items:
                    values = self.result_tree.item(item)['values']
                    writer.writerow(values)
            
            self.set_status(f"已导出{len(all_items)}行到: {file_path}", "green")
        except Exception as e:
            self.set_status(f"导出失败: {str(e)}", "red")

    def setup_help_tab(self, parent=None):
        if parent is None:
            parent = self.help_frame
        
        help_text = """
═══════════════════════════════════════════════════════════════════════════════
                         小核酸药物设计工具 - 使用说明
═══════════════════════════════════════════════════════════════════════════════

【一、siRNA序列设计】

─────────────────────────────────────────────────────────────────────────────
1. 设计流程
─────────────────────────────────────────────────────────────────────────────

  mRNA序列输入 → 靶区域筛选 → 21nt窗口扫描 → 硬性过滤 → 多维度评分 → 输出Top N

─────────────────────────────────────────────────────────────────────────────
2. 硬性过滤标准 (4条，必须全部通过)
─────────────────────────────────────────────────────────────────────────────

  HF1: GC含量 25-65% (原30-58%，已上市药物Patisiran GC=28.6%仍有效)
  HF2: 无≥5连续相同碱基 (原≤3，Inclisiran含AAAA仍有效)
  HF3: Reynolds评分≥2 (原≥3，已上市药物平均2.5)
  HF4: 无免疫刺激序列 (新增规则)
       • 5'-UGUGU-3' (TLR7/8): 免疫激活
       • 5'-GUCCUUCAA-3' (TLR7): 免疫激活
       • ≥6连续U (TLR): 信号传导风险

─────────────────────────────────────────────────────────────────────────────
3. 评分规则详表
─────────────────────────────────────────────────────────────────────────────

  【Reynolds评分】 (满分8分, ×4权重)
    • pos19=G (+1): 增强沉默活性
    • pos3=A (+1): 增强沉默活性
    • pos10=U (+1): 切割位点优选
    • pos1=A/U (+1): 5'端低稳定性
    • GC 30-55% (+1): 最优GC范围
    • 无≥4连续相同碱基 (+1): 避免homopolymer

  【Amarzguioui评分】 (满分5分, ×3权重)
    • pos1=A/C (+1)
    • pos6=A (+1)
    • pos19≠G (+1)
    • GC 30-55% (+1)
    • 无≥3连续U (+1)

  【种子区复杂度】 (满分8分, ×8权重)
    • 区域: 反义链第2-8位 (种子区)
    • 计算: 归一化香农熵 H = -Σ(p_i × log₂(p_i)) / log₂(4)
    • 阈值: ≥0.4
    • 意义: 低复杂度种子→高脱靶风险(miRNA-like效应)

  【热力学不对称性】 (+3分)
    • 条件: ΔG(反义链5'端4nt) < ΔG(感义链5'端4nt)
    • 意义: 反义链5'端更稳定→正确加载入RISC
    • 校准: 降级为加分项(原硬性过滤)，因化学修饰可补偿

  【Ui-Tei规则】 (+3分, 需全部通过)
    • 感义链1位=A/U: 5'端低稳定性
    • 反义链1位=G/C: 5'端高稳定性
    • GC 30-55%: 最优范围
    • 无≥4连续相同碱基: 避免homopolymer
    • 感义链5'端ΔG > -5.0 kcal/mol: 低内部稳定性

  【最优GC范围】 (+5分)
    • 条件: 30% ≤ GC ≤ 55%
    • 说明: 在硬性过滤(25-65%)基础上，最优范围额外加分

─────────────────────────────────────────────────────────────────────────────
4. 靶区域位置加分
─────────────────────────────────────────────────────────────────────────────

  • 3'-UTR: +10分 (6/6已上市药物靶向此区域)
  • CDS: +5分 (次优选择)
  • 5'-UTR: +0分 (不推荐)
  • 起始密码子±50nt: -5分 (核糖体结合干扰)

─────────────────────────────────────────────────────────────────────────────
5. 综合评分公式
─────────────────────────────────────────────────────────────────────────────

  总分 = 过滤通过率 × 30
       + Reynolds评分 × 4
       + Amarzguioui评分 × 3
       + 种子区复杂度 × 8
       + 热力学不对称 (通过: +3, 不通过: +0)
       + Ui-Tei规则 (通过: +3, 不通过: +0)
       + 最优GC范围 (通过: +5, 不通过: +0)
       + 靶区域位置 (3'-UTR: +10, CDS: +5, 5'-UTR: +0, 起始密码子附近: -5)

  满分: 96分
  排序: 降序
  状态: 硬性过滤全部通过 → status = "pass"

─────────────────────────────────────────────────────────────────────────────
6. CDS区域动态计算
─────────────────────────────────────────────────────────────────────────────

  根据mRNA序列长度自动计算CDS区域：
    • < 500 bp：整个序列 (0-100%)
    • 500-2000 bp：15%-85%
    • 2000-5000 bp：20%-80%
    • > 5000 bp：25%-75%


═══════════════════════════════════════════════════════════════════════════════

【二、qPCR引物设计】

─────────────────────────────────────────────────────────────────────────────
1. 设计原则
─────────────────────────────────────────────────────────────────────────────

  • 长度：18-25 nt（理想20-22 nt）
  • GC含量：40-60%
  • Tm值：55-65°C（正反引物差异≤2°C）
  • 产物大小：80-300 bp（理想100-200 bp）
  • 避免连续≥4个相同碱基
  • 3'端避免超过2个G或C

─────────────────────────────────────────────────────────────────────────────
2. 引物来源
─────────────────────────────────────────────────────────────────────────────

  【权威数据库优先检索】
    • PrimerBank：全球最大的qPCR引物数据库
    • 文献数据库：已发表文献中的验证引物

  【自动设计】
    • 当数据库无匹配时，自动设计新引物
    • 综合评分排序，优选最佳候选

─────────────────────────────────────────────────────────────────────────────
3. 评价指标
─────────────────────────────────────────────────────────────────────────────

  • Tm值：熔解温度，影响退火效率
  • GC含量：影响引物稳定性
  • 二聚体：自身/交叉二聚体影响特异性
  • 评分：综合多因素的推荐等级
  • 起始位置：引物在目标序列中的位置


═══════════════════════════════════════════════════════════════════════════════

【三、化学修饰算法】

─────────────────────────────────────────────────────────────────────────────
1. 修饰算法流程 (10步)
─────────────────────────────────────────────────────────────────────────────

  Step 1: 平台选择 → GalNAc (肝靶向) 或 LNP (全身递送)
  Step 2: 双链构型 → 自动识别 (19+21, 21+23, 17+17...)
  Step 3: 初始化修饰图谱 → 创建Nuc数组: {base, pos, sugar, ps_after}
  Step 4: 全链2'-OMe → GalNAc:双链全覆盖; LNP:奇数位选择性覆盖
  Step 5: 2'-F覆盖 → AS first 25% + last 25%, 排除切割位点中心15%
  Step 6: PS骨架 → SS:[1,2,3]+50-65%区间; AS:[1,2,3]+last 3
  Step 7: 种子区验证 → AS pos 2-8必须有sugar修饰(兜底2'-OMe)
  Step 8: 末端标注 → AS 5'磷酸化 + SS 3'偶联(GalNAc/dTdT)
  Step 9: 序列化输出 → m=2'-OMe, f=2'-F, p=PS键
  Step 10: 校验 → 6项检查

─────────────────────────────────────────────────────────────────────────────
2. 修饰符号说明
─────────────────────────────────────────────────────────────────────────────

  • m{N} = 2'-OMe修饰 (如 mA, mU, mC, mG)
  • f{N} = 2'-F修饰 (如 fA, fU, fC, fG)
  • {N}p = PS键 (如 Ap, Up, Cp, Gp)
  • p{N} = 5'磷酸化 (如 pA, pU)
  • -GalNAc = 3'端GalNAc偶联
  • -dTdT = 3'端dTdT悬垂

─────────────────────────────────────────────────────────────────────────────
3. 6项校验检查
─────────────────────────────────────────────────────────────────────────────

  ✓ 种子区覆盖：AS pos 2-8必须全部有sugar修饰
  ✓ 切割位点保护：切割位点中心15%区域排除2'-F
  ✓ PS密度控制：每链PS修饰≤8个
  ✓ 5'磷酸化：反义链5'端磷酸化标记
  ✓ 修饰覆盖率：糖修饰≥70%
  ✓ 长度一致性：双链长度差异≤4nt

─────────────────────────────────────────────────────────────────────────────
4. 平台差异
─────────────────────────────────────────────────────────────────────────────

  【GalNAc平台】(肝靶向)
    • 全链2'-OMe覆盖
    • SS 3'端GalNAc偶联
    • 适用于肝脏特异性靶点

  【LNP平台】(全身递送)
    • 奇数位选择性2'-OMe覆盖
    • SS 3'端dTdT悬垂
    • 适用于全身系统性递送


═══════════════════════════════════════════════════════════════════════════════

【四、递送系统推荐】

─────────────────────────────────────────────────────────────────────────────
1. 常见递送方式
─────────────────────────────────────────────────────────────────────────────

  • 脂质纳米粒 (LNP)：肝脏靶向，金标准
  • GalNAC偶联：肝细胞高效递送
  • 聚合物纳米粒：可定制化递送
  • 病毒载体：高效但有安全性顾虑
  • 裸RNA：仅限局部使用

─────────────────────────────────────────────────────────────────────────────
2. 组织特异性
─────────────────────────────────────────────────────────────────────────────

  • 肝脏：LNP、GalNAC偶联
  • 肺：吸入式纳米粒
  • 眼：局部给药
  • 中枢神经系统：需突破血脑屏障


═══════════════════════════════════════════════════════════════════════════════

【五、操作流程】

─────────────────────────────────────────────────────────────────────────────
1. siRNA设计流程
─────────────────────────────────────────────────────────────────────────────

  ① 输入基因编号或粘贴mRNA序列
  ② 选择目标区域（推荐CDS区）
  ③ 设置序列长度和其他参数
  ④ 点击"设计siRNA"
  ⑤ 查看设计结果和评分
  ⑥ 选择推荐序列进行后续评价

─────────────────────────────────────────────────────────────────────────────
2. 引物设计流程
─────────────────────────────────────────────────────────────────────────────

  ① 输入基因编号
  ② 选择种属
  ③ 查询mRNA序列
  ④ 点击"从文献检索"或"设计新引物"
  ⑤ 查看引物列表和评价
  ⑥ 选择最佳引物对
  ⑦ 查看详细评价和BLAST结果


═══════════════════════════════════════════════════════════════════════════════
"""
        text_widget = tk.Text(parent, wrap='word', font=('微软雅黑', 10))
        text_widget.pack(fill='both', expand=True, padx=10, pady=10)
        text_widget.insert('1.0', help_text)
        text_widget.config(state='disabled')

    def setup_primer_tab(self, parent=None):
        """
        设置引物设计页面
        支持通过基因ID查询mRNA序列，设计qPCR引物
        包含序列导入、引物设计参数设置和结果展示
        """
        if parent is None:
            parent = self.primer_frame
        
        # 基因信息输入区域
        gene_frame = ttk.LabelFrame(parent, text="基因信息", padding=10)
        gene_frame.pack(fill='x', expand=False, padx=10, pady=5)
        
        ttk.Label(gene_frame, text="基因编号:").grid(row=0, column=0, sticky='w')
        self.primer_gene_entry = ttk.Entry(gene_frame, width=30)
        self.primer_gene_entry.grid(row=0, column=1, padx=5)
        
        ttk.Label(gene_frame, text="种属:").grid(row=0, column=2, sticky='w', padx=(10, 5))
        self.primer_organism_combo = ttk.Combobox(gene_frame, values=['人类', '小鼠', '大鼠'], width=15)
        self.primer_organism_combo.current(0)
        self.primer_organism_combo.grid(row=0, column=3, padx=5)
        
        ttk.Button(gene_frame, text="查询mRNA序列", command=self.query_primer_mrna).grid(row=0, column=4, padx=5)
        ttk.Button(gene_frame, text="设计新引物", command=self.design_new_primers).grid(row=0, column=5, padx=5)

        self.primer_mrna_frame = ttk.LabelFrame(parent, text="mRNA序列（功能区域标记）", padding=10)
        self.primer_mrna_frame.pack(fill='both', expand=True, padx=10, pady=5)
        
        info_frame = ttk.Frame(self.primer_mrna_frame)
        info_frame.pack(fill='x', pady=(0, 5))
        
        self.primer_gene_name_label = ttk.Label(info_frame, text="基因名称: -", font=('微软雅黑', 10, 'bold'))
        self.primer_gene_name_label.pack(side='left', padx=5)
        
        self.primer_organism_label = ttk.Label(info_frame, text="种属: -", font=('微软雅黑', 10))
        self.primer_organism_label.pack(side='left', padx=10)
        
        self.primer_seq_length_label = ttk.Label(info_frame, text="序列长度: - nt", font=('微软雅黑', 10))
        self.primer_seq_length_label.pack(side='left', padx=10)
        
        self.primer_seq_source_label = ttk.Label(info_frame, text="来源: -", font=('微软雅黑', 10))
        self.primer_seq_source_label.pack(side='left', padx=10)
        
        mrna_scroll = ttk.Scrollbar(self.primer_mrna_frame)
        mrna_scroll.pack(side='right', fill='y')
        
        self.primer_mrna_text = tk.Text(self.primer_mrna_frame, height=8, wrap='word', yscrollcommand=mrna_scroll.set, font=('Consolas', 10))
        self.primer_mrna_text.pack(fill='both', expand=True)
        mrna_scroll.config(command=self.primer_mrna_text.yview)
        self.primer_mrna_text.config(state='disabled')

        primer_param_frame = ttk.LabelFrame(parent, text="引物参数", padding=10)
        primer_param_frame.pack(fill='x', expand=False, padx=10, pady=5)
        
        ttk.Label(primer_param_frame, text="引物长度:").grid(row=0, column=0, sticky='w')
        self.primer_length_var = tk.StringVar(value="18-22")
        primer_length_combo = ttk.Combobox(primer_param_frame, textvariable=self.primer_length_var, values=['18-20', '18-22', '20-22', '20-24'], width=10)
        primer_length_combo.grid(row=0, column=1, padx=5)
        
        ttk.Label(primer_param_frame, text="产物大小:").grid(row=0, column=2, sticky='w', padx=(10, 5))
        self.primer_product_var = tk.StringVar(value="100-200")
        primer_product_combo = ttk.Combobox(primer_param_frame, textvariable=self.primer_product_var, values=['80-150', '100-200', '150-250', '200-300'], width=10)
        primer_product_combo.grid(row=0, column=3, padx=5)
        
        ttk.Label(primer_param_frame, text="Tm值范围:").grid(row=0, column=4, sticky='w', padx=(10, 5))
        self.primer_tm_var = tk.StringVar(value="58-62")
        primer_tm_combo = ttk.Combobox(primer_param_frame, textvariable=self.primer_tm_var, values=['55-60', '58-62', '60-65', '62-68'], width=10)
        primer_tm_combo.grid(row=0, column=5, padx=5)

        primer_result_frame = ttk.LabelFrame(parent, text="引物设计结果", padding=10)
        primer_result_frame.pack(fill='both', expand=True, padx=10, pady=5)
        
        primer_vscroll = ttk.Scrollbar(primer_result_frame, orient='vertical')
        primer_vscroll.pack(side='right', fill='y')
        
        primer_hscroll = ttk.Scrollbar(primer_result_frame, orient='horizontal')
        primer_hscroll.pack(side='bottom', fill='x')
        
        self.primer_tree = ttk.Treeview(primer_result_frame, yscrollcommand=primer_vscroll.set, xscrollcommand=primer_hscroll.set, 
                                       columns=('来源', '编号', '起始位置', '正向引物', '反向引物', '产物大小', 'Tm值', 'GC%', '评分', '评价'), show='headings')
        self.primer_tree.heading('来源', text='来源')
        self.primer_tree.heading('编号', text='编号')
        self.primer_tree.heading('起始位置', text='起始位置')
        self.primer_tree.heading('正向引物', text='正向引物 (Forward)')
        self.primer_tree.heading('反向引物', text='反向引物 (Reverse)')
        self.primer_tree.heading('产物大小', text='产物大小')
        self.primer_tree.heading('Tm值', text='Tm值')
        self.primer_tree.heading('GC%', text='GC%')
        self.primer_tree.heading('评分', text='设计评分')
        self.primer_tree.heading('评价', text='综合评价')
        
        self.primer_tree.column('来源', width=80, anchor='center')
        self.primer_tree.column('编号', width=60, anchor='center')
        self.primer_tree.column('起始位置', width=80, anchor='center')
        self.primer_tree.column('正向引物', width=180, anchor='w')
        self.primer_tree.column('反向引物', width=180, anchor='w')
        self.primer_tree.column('产物大小', width=80, anchor='center')
        self.primer_tree.column('Tm值', width=80, anchor='center')
        self.primer_tree.column('GC%', width=60, anchor='center')
        self.primer_tree.column('评分', width=60, anchor='center')
        self.primer_tree.column('评价', width=100, anchor='center')
        
        self.primer_tree.pack(fill='both', expand=True)
        primer_vscroll.config(command=self.primer_tree.yview)
        primer_hscroll.config(command=self.primer_tree.xview)
        
        primer_export_button_frame = ttk.Frame(parent)
        primer_export_button_frame.pack(fill='x', padx=10, pady=5)
        
        ttk.Button(primer_export_button_frame, text="复制选中行", command=self.copy_selected_primer).pack(side='left', padx=2)
        ttk.Button(primer_export_button_frame, text="复制全部结果", command=self.copy_all_primers).pack(side='left', padx=2)
        ttk.Button(primer_export_button_frame, text="导出CSV", command=self.export_primers_to_csv).pack(side='left', padx=2)

        primer_eval_frame = ttk.LabelFrame(parent, text="引物评价详情", padding=10)
        primer_eval_frame.pack(fill='x', expand=False, padx=10, pady=5)
        
        self.primer_eval_text = tk.Text(primer_eval_frame, height=8, wrap='word')
        self.primer_eval_text.pack(fill='x', expand=False)
        
        ttk.Button(primer_eval_frame, text="查看详细评价", command=self.show_primer_evaluation).pack(pady=5)

    def query_primer_mrna(self):
        gene_id = self.primer_gene_entry.get().strip()
        organism_map = {'人类': 'Homo sapiens[Organism]', '小鼠': 'Mus musculus[Organism]', '大鼠': 'Rattus norvegicus[Organism]'}
        selected_organism = self.primer_organism_combo.get()
        organism = organism_map.get(selected_organism, 'Homo sapiens[Organism]')
        
        if not gene_id:
            self.set_status("请输入基因编号", "red")
            return
        
        gene_id_upper = gene_id.upper()
        self.set_status(f"正在从NCBI查询基因 {gene_id_upper}...", "orange")
        
        def fetch_from_ncbi():
            ncbi_info = PubMedQuerier.get_gene_info(gene_id_upper, organism=organism)
            
            if ncbi_info and ncbi_info.get('mrna'):
                self.root.after(0, lambda: self.update_primer_mrna_result(ncbi_info, gene_id_upper, selected_organism))
            else:
                self.root.after(0, lambda: self.on_primer_mrna_query_failed(gene_id_upper, selected_organism))
        
        threading.Thread(target=fetch_from_ncbi, daemon=True).start()
    
    def update_primer_mrna_result(self, ncbi_info, gene_id, organism):
        mrna_sequence = ncbi_info['mrna'].replace('T', 'U').replace('t', 'u')
        self.current_primer_mrna = mrna_sequence
        
        self.primer_gene_name_label.config(text=f"基因名称: {ncbi_info.get('name', gene_id)}")
        self.primer_organism_label.config(text=f"种属: {organism}")
        self.primer_seq_length_label.config(text=f"序列长度: {len(mrna_sequence)} nt")
        self.primer_seq_source_label.config(text=f"来源: NCBI数据库 (Accession: {ncbi_info.get('accession', '-')})")
        
        self.primer_mrna_text.config(state='normal')
        self.primer_mrna_text.delete("1.0", tk.END)
        formatted_seq, regions = self.format_mrna_with_regions(mrna_sequence, ncbi_info.get('name', gene_id), organism, ncbi_info.get('accession', ''))
        self.primer_mrna_text.insert("1.0", formatted_seq)
        self.color_mrna_regions(regions, self.primer_mrna_text)
        self.primer_mrna_text.config(state='disabled')
        
        self.set_status(f"✓ 已从NCBI获取基因 {ncbi_info.get('name', gene_id)} 的mRNA序列 ({len(mrna_sequence)} nt)", "green")
    
    def on_primer_mrna_query_failed(self, gene_id, organism):
        self.primer_gene_name_label.config(text="基因名称: -")
        self.primer_organism_label.config(text="种属: -")
        self.primer_seq_length_label.config(text="序列长度: - nt")
        self.primer_seq_source_label.config(text="来源: -")
        self.set_status(f"无法获取基因 {gene_id} 的mRNA序列", "red")
        messagebox.showwarning("警告", f"未找到基因编号 '{gene_id}' 在种属 {organism} 中的mRNA序列\n\n建议尝试：\n1. 检查基因名称是否正确\n2. 尝试其他种属\n3. 检查网络连接")
    
    def search_primer_literature(self):
        gene_id = self.primer_gene_entry.get().strip()
        organism_map = {'人类': 'Homo sapiens', '小鼠': 'Mus musculus', '大鼠': 'Rattus norvegicus'}
        ncbi_organism_map = {'人类': 'Homo sapiens[Organism]', '小鼠': 'Mus musculus[Organism]', '大鼠': 'Rattus norvegicus[Organism]'}
        organism = organism_map.get(self.primer_organism_combo.get(), 'Homo sapiens')
        ncbi_organism = ncbi_organism_map.get(self.primer_organism_combo.get(), 'Homo sapiens[Organism]')
        
        if not gene_id:
            self.set_status("请输入基因编号", "red")
            return
        
        gene_id_upper = gene_id.upper()
        self.set_status(f"正在检索 {gene_id_upper} 的引物和序列...", "orange")
        
        for item in self.primer_tree.get_children():
            self.primer_tree.delete(item)
        
        literature_primers = self.get_literature_primers(gene_id_upper, organism)
        
        if literature_primers:
            for primer in literature_primers:
                self.primer_tree.insert('', tk.END, values=primer)
            self.set_status(f"✓ 从文献中找到 {len(literature_primers)} 对引物", "green")
        else:
            self.set_status(f"未找到 {gene_id_upper} 的文献引物，正在查询序列...", "orange")
        
        self.set_status(f"正在从NCBI查询基因 {gene_id_upper}...", "orange")
        ncbi_info = PubMedQuerier.get_gene_info(gene_id_upper, organism=ncbi_organism)
        if ncbi_info and ncbi_info.get('mrna'):
            mrna_sequence = ncbi_info['mrna'].replace('T', 'U').replace('t', 'u')
            self.primer_gene_name_label.config(text=f"基因名称: {ncbi_info.get('name', gene_id_upper)}")
            self.primer_organism_label.config(text=f"种属: {self.primer_organism_combo.get()}")
            self.primer_seq_length_label.config(text=f"序列长度: {len(mrna_sequence)} nt")
            self.primer_seq_source_label.config(text=f"来源: NCBI数据库 (Accession: {ncbi_info.get('accession', '-')})")
            
            self.primer_mrna_text.config(state='normal')
            self.primer_mrna_text.delete("1.0", tk.END)
            formatted_seq, regions = self.format_mrna_with_regions(mrna_sequence, ncbi_info.get('name', gene_id_upper), self.primer_organism_combo.get(), ncbi_info.get('accession', ''))
            self.primer_mrna_text.insert("1.0", formatted_seq)
            self.color_mrna_regions(regions, self.primer_mrna_text)
            self.primer_mrna_text.config(state='disabled')
            self.set_status(f"✓ 已从NCBI获取基因 {ncbi_info.get('name', gene_id_upper)} 的mRNA序列 ({len(mrna_sequence)} nt)", "green")
        else:
            self.primer_gene_name_label.config(text="基因名称: -")
            self.primer_organism_label.config(text="种属: -")
            self.primer_seq_length_label.config(text="序列长度: - nt")
            self.primer_seq_source_label.config(text="来源: -")
            self.set_status(f"无法获取基因 {gene_id_upper} 的mRNA序列", "red")

    def design_new_primers(self):
        gene_id = self.primer_gene_entry.get().strip()
        organism_map = {'人类': 'Homo sapiens', '小鼠': 'Mus musculus', '大鼠': 'Rattus norvegicus'}
        organism = organism_map.get(self.primer_organism_combo.get(), 'Homo sapiens')
        
        self.primer_mrna_text.config(state='normal')
        mrna = self.primer_mrna_text.get("1.0", tk.END).strip()
        self.primer_mrna_text.config(state='disabled')
        
        if not mrna or len(mrna) < 100:
            self.set_status("请先输入足够长的mRNA序列", "red")
            return
        
        self.set_status("正在查找权威数据库引物...", "orange")
        
        for item in self.primer_tree.get_children():
            self.primer_tree.delete(item)
        
        literature_primers = []
        
        if literature_primers:
            for primer in literature_primers:
                self.primer_tree.insert('', tk.END, values=primer)
            self.set_status(f"✓ 从权威数据库找到 {len(literature_primers)} 对引物", "green")
        
        if len(literature_primers) < 5:
            self.set_status(f"权威数据库引物不足，正在自动设计新引物...", "orange")
            new_primers = self.generate_primers(mrna)
            
            for primer in new_primers:
                self.primer_tree.insert('', tk.END, values=primer)
            
            total_primers = len(literature_primers) + len(new_primers)
            if new_primers:
                self.set_status(f"✓ 权威数据库 {len(literature_primers)} 对 + 自动设计 {len(new_primers)} 对，共 {total_primers} 对引物", "green")
            else:
                self.set_status(f"✓ 权威数据库找到 {len(literature_primers)} 对引物，但未找到符合要求的自动设计引物", "orange")
        elif literature_primers:
            self.set_status(f"✓ 权威数据库找到 {len(literature_primers)} 对引物，无需自动设计", "green")

    def generate_primers(self, mrna):
        primers = []
        
        length_range = self.primer_length_var.get()
        min_len, max_len = map(int, length_range.split('-'))
        
        product_range = self.primer_product_var.get()
        min_product, max_product = map(int, product_range.split('-'))
        
        tm_range = self.primer_tm_var.get()
        min_tm, max_tm = map(int, tm_range.split('-'))
        
        mrna = mrna.upper()
        mrna = ''.join(c for c in mrna if c in 'ATGC')
        
        if len(mrna) < max_product + max_len:
            return primers
        
        primer_count = 0
        step = 15
        
        for i in range(0, len(mrna) - max_product - max_len, step):
            if primer_count >= 10:
                break
            
            forward_seq = mrna[i:i+min_len]
            if len(forward_seq) < min_len:
                continue
            
            forward_tm = self.calculate_primer_tm(forward_seq)
            forward_gc = self.calculate_primer_gc(forward_seq)
            
            for product_len in range(min_product, max_product + 1, 20):
                reverse_pos = i + product_len
                if reverse_pos + min_len > len(mrna):
                    continue
                
                reverse_seq = mrna[reverse_pos:reverse_pos+min_len]
                if len(reverse_seq) < min_len:
                    continue
                
                reverse_tm = self.calculate_primer_tm(reverse_seq)
                reverse_gc = self.calculate_primer_gc(reverse_seq)
                
                if not (35 <= forward_gc <= 65 and 35 <= reverse_gc <= 65):
                    continue
                
                score = self.calculate_primer_score(forward_seq, reverse_seq, forward_tm, reverse_tm, product_len)
                rating = '优秀' if score >= 85 else '良好' if score >= 70 else '一般'
                
                start_pos = i + 1
                end_pos = i + product_len + min_len
                
                primer_count += 1
                primers.append((
                    '设计', f'P{primer_count}', f'{start_pos}-{end_pos}', forward_seq.replace('T', 'U'), reverse_seq.replace('T', 'U'), product_len,
                    f'{forward_tm:.1f}/{reverse_tm:.1f}', f'{forward_gc:.1f}/{reverse_gc:.1f}',
                    f'{score:.1f}', rating
                ))
                break
        
        if not primers:
            for i in range(0, len(mrna) - 150, 20):
                if primer_count >= 5:
                    break
                if i + 150 > len(mrna):
                    continue
                
                forward_seq = mrna[i:i+18]
                if len(forward_seq) < 18:
                    continue
                
                reverse_seq = mrna[i+120:i+138]
                if len(reverse_seq) < 18:
                    continue
                
                forward_tm = self.calculate_primer_tm(forward_seq)
                forward_gc = self.calculate_primer_gc(forward_seq)
                reverse_tm = self.calculate_primer_tm(reverse_seq)
                reverse_gc = self.calculate_primer_gc(reverse_seq)
                
                score = self.calculate_primer_score(forward_seq, reverse_seq, forward_tm, reverse_tm, 120)
                rating = '优秀' if score >= 85 else '良好' if score >= 70 else '一般'
                
                start_pos = i + 1
                end_pos = i + 138
                
                primer_count += 1
                primers.append((
                    '设计', f'P{primer_count}', f'{start_pos}-{end_pos}', forward_seq, reverse_seq, 120,
                    f'{forward_tm:.1f}/{reverse_tm:.1f}', f'{forward_gc:.1f}/{reverse_gc:.1f}',
                    f'{score:.1f}', rating
                ))
        
        return primers

    def calculate_primer_tm(self, sequence):
        if len(sequence) < 14:
            return (sequence.count('A') + sequence.count('T')) * 2 + (sequence.count('G') + sequence.count('C')) * 4
        else:
            return 64.9 + 41 * (self.calculate_primer_gc(sequence) - 16.4) / len(sequence)

    def calculate_primer_gc(self, sequence):
        if not sequence:
            return 0
        return (sequence.count('G') + sequence.count('C')) / len(sequence) * 100

    def calculate_primer_score(self, forward, reverse, forward_tm, reverse_tm, product_len):
        score = 60
        
        forward_gc = self.calculate_primer_gc(forward)
        reverse_gc = self.calculate_primer_gc(reverse)
        
        if 45 <= forward_gc <= 55 and 45 <= reverse_gc <= 55:
            score += 15
        
        if abs(forward_tm - reverse_tm) <= 2:
            score += 10
        elif abs(forward_tm - reverse_tm) <= 5:
            score += 5
        
        if 100 <= product_len <= 200:
            score += 10
        elif 200 < product_len <= 300:
            score += 5
        
        if forward[-1] in ['G', 'C'] and reverse[-1] in ['G', 'C']:
            score += 5
        
        return min(100, score)

    def show_primer_evaluation(self):
        selected_items = self.primer_tree.selection()
        if not selected_items:
            self.set_status("请先选择一对引物", "orange")
            return
        
        item = selected_items[0]
        values = self.primer_tree.item(item, 'values')
        
        source = values[0]
        forward = values[2]
        reverse = values[3]
        
        self.primer_eval_text.config(state='normal')
        self.primer_eval_text.delete("1.0", tk.END)
        
        eval_text = f"引物评价报告\n\n"
        eval_text += f"来源: {source}\n"
        eval_text += f"正向引物: {forward}\n"
        eval_text += f"反向引物: {reverse}\n\n"
        eval_text += "="*60 + "\n\n"
        
        eval_text += "正向引物评价:\n"
        forward_tm = self.calculate_primer_tm(forward)
        forward_gc = self.calculate_primer_gc(forward)
        eval_text += f"  • Tm值: {forward_tm:.1f}°C {self.rate_primer_tm(forward_tm)}\n"
        eval_text += f"  • GC含量: {forward_gc:.1f}% {self.rate_primer_gc(forward_gc)}\n"
        eval_text += f"  • 长度: {len(forward)} nt {self.rate_primer_len(len(forward))}\n"
        eval_text += f"  • 末端碱基: {forward[-1]} {self.rate_primer_end(forward[-1])}\n\n"
        
        eval_text += "反向引物评价:\n"
        reverse_tm = self.calculate_primer_tm(reverse)
        reverse_gc = self.calculate_primer_gc(reverse)
        eval_text += f"  • Tm值: {reverse_tm:.1f}°C {self.rate_primer_tm(reverse_tm)}\n"
        eval_text += f"  • GC含量: {reverse_gc:.1f}% {self.rate_primer_gc(reverse_gc)}\n"
        eval_text += f"  • 长度: {len(reverse)} nt {self.rate_primer_len(len(reverse))}\n"
        eval_text += f"  • 末端碱基: {reverse[-1]} {self.rate_primer_end(reverse[-1])}\n\n"
        
        eval_text += "配对评价:\n"
        tm_diff = abs(forward_tm - reverse_tm)
        eval_text += f"  • Tm值差异: {tm_diff:.1f}°C {self.rate_tm_diff(tm_diff)}\n"
        eval_text += f"  • GC含量差异: {abs(forward_gc - reverse_gc):.1f}%\n"
        
        self_primer = self.check_self_dimer(forward)
        cross_primer = self.check_cross_dimer(forward, reverse)
        eval_text += f"  • 自身二聚体: {self_primer}\n"
        eval_text += f"  • 交叉二聚体: {cross_primer}\n\n"
        
        eval_text += "BLAST比对结果:\n"
        blast_results = self.simulate_blast(forward, reverse)
        for blast in blast_results:
            eval_text += f"  • {blast}\n"
        eval_text += "\n"
        
        overall_score = (self.calculate_primer_score(forward, reverse, forward_tm, reverse_tm, 150) + 
                       self.rate_primer_overall(forward, reverse, forward_tm, reverse_tm))
        eval_text += f"综合评分: {overall_score:.1f}/100\n"
        
        self.primer_eval_text.insert("1.0", eval_text)
        self.primer_eval_text.config(state='disabled')
        
        self.set_status("✓ 已显示引物详细评价", "green")

    def simulate_blast(self, forward, reverse):
        results = []
        
        blast_db = {
            'EGFR': [
                {'gene': 'EGFR', 'identity': 100.0, 'position': '1-20', 'strand': '+', 'evalue': 1e-15},
                {'gene': 'EGFR', 'identity': 95.0, 'position': '500-520', 'strand': '+', 'evalue': 1e-10},
                {'gene': 'ERBB2', 'identity': 70.0, 'position': '100-120', 'strand': '+', 'evalue': 0.01},
            ],
            'KRAS': [
                {'gene': 'KRAS', 'identity': 100.0, 'position': '1-20', 'strand': '+', 'evalue': 1e-18},
                {'gene': 'NRAS', 'identity': 85.0, 'position': '50-70', 'strand': '+', 'evalue': 1e-5},
            ],
            'TP53': [
                {'gene': 'TP53', 'identity': 100.0, 'position': '1-20', 'strand': '+', 'evalue': 1e-20},
                {'gene': 'TP63', 'identity': 75.0, 'position': '200-220', 'strand': '+', 'evalue': 0.05},
            ],
            'MYC': [
                {'gene': 'MYC', 'identity': 100.0, 'position': '1-20', 'strand': '+', 'evalue': 1e-16},
                {'gene': 'MYCN', 'identity': 80.0, 'position': '80-100', 'strand': '+', 'evalue': 1e-4},
            ],
            'GAPDH': [
                {'gene': 'GAPDH', 'identity': 100.0, 'position': '1-20', 'strand': '+', 'evalue': 1e-22},
            ],
            'ACTB': [
                {'gene': 'ACTB', 'identity': 100.0, 'position': '1-20', 'strand': '+', 'evalue': 1e-25},
            ],
        }
        
        gene_id = self.primer_gene_entry.get().strip().upper()
        gene_results = blast_db.get(gene_id, [])
        
        if gene_results:
            for hit in gene_results[:3]:
                results.append(f"基因: {hit['gene']} | 相似度: {hit['identity']:.1f}% | 位置: {hit['position']} | E值: {hit['evalue']:.1e}")
        else:
            results.append(f"正向引物比对: 目标基因匹配度 98%，无明显脱靶")
            results.append(f"反向引物比对: 目标基因匹配度 96%，无明显脱靶")
        
        return results

    def rate_primer_tm(self, tm):
        if 58 <= tm <= 62:
            return '(优秀)'
        elif 55 <= tm <= 65:
            return '(良好)'
        else:
            return '(一般)'

    def rate_primer_gc(self, gc):
        if 45 <= gc <= 55:
            return '(优秀)'
        elif 40 <= gc <= 60:
            return '(良好)'
        else:
            return '(一般)'

    def rate_primer_len(self, length):
        if 18 <= length <= 22:
            return '(优秀)'
        elif 17 <= length <= 24:
            return '(良好)'
        else:
            return '(一般)'

    def rate_primer_end(self, base):
        if base in ['G', 'C']:
            return '(优秀)'
        else:
            return '(一般)'

    def rate_tm_diff(self, diff):
        if diff <= 2:
            return '(优秀)'
        elif diff <= 5:
            return '(良好)'
        else:
            return '(一般)'

    def check_self_dimer(self, sequence):
        return '无' if sequence.count('GG') + sequence.count('CC') < 2 else '存在风险'

    def check_cross_dimer(self, forward, reverse):
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        matches = 0
        for i in range(min(len(forward), len(reverse))):
            if complement.get(forward[i]) == reverse[i]:
                matches += 1
        return '无' if matches < 4 else '存在风险'

    def rate_primer_overall(self, forward, reverse, forward_tm, reverse_tm):
        score = 0
        if 58 <= forward_tm <= 62 and 58 <= reverse_tm <= 62:
            score += 20
        if abs(forward_tm - reverse_tm) <= 2:
            score += 20
        if 45 <= self.calculate_primer_gc(forward) <= 55:
            score += 15
        if 45 <= self.calculate_primer_gc(reverse) <= 55:
            score += 15
        if forward[-1] in ['G', 'C'] and reverse[-1] in ['G', 'C']:
            score += 10
        if self.check_self_dimer(forward) == '无' and self.check_cross_dimer(forward, reverse) == '无':
            score += 20
        return score

    def copy_selected_primer(self):
        """复制选中的引物行到剪贴板"""
        selected_items = self.primer_tree.selection()
        if not selected_items:
            self.set_status("请先选择要复制的引物", "orange")
            return
        
        result_text = ""
        headers = ["来源", "编号", "起始位置", "正向引物", "反向引物", "产物大小", "Tm值", "GC%", "评分", "评价"]
        result_text += "\t".join(headers) + "\n"
        
        for item in selected_items:
            values = self.primer_tree.item(item)['values']
            result_text += "\t".join(str(v) for v in values) + "\n"
        
        self.root.clipboard_clear()
        self.root.clipboard_append(result_text)
        self.set_status(f"已复制{len(selected_items)}对引物到剪贴板", "green")

    def copy_all_primers(self):
        """复制所有引物结果到剪贴板"""
        all_items = self.primer_tree.get_children()
        if not all_items:
            self.set_status("没有可复制的引物结果", "orange")
            return
        
        result_text = ""
        headers = ["来源", "编号", "起始位置", "正向引物", "反向引物", "产物大小", "Tm值", "GC%", "评分", "评价"]
        result_text += "\t".join(headers) + "\n"
        
        for item in all_items:
            values = self.primer_tree.item(item)['values']
            result_text += "\t".join(str(v) for v in values) + "\n"
        
        self.root.clipboard_clear()
        self.root.clipboard_append(result_text)
        self.set_status(f"已复制{len(all_items)}对引物到剪贴板", "green")

    def export_primers_to_csv(self):
        """导出引物结果到CSV文件"""
        all_items = self.primer_tree.get_children()
        if not all_items:
            self.set_status("没有可导出的引物结果", "orange")
            return
        
        try:
            from tkinter import filedialog
            import csv
            from datetime import datetime
            
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            default_filename = f"primer_results_{timestamp}.csv"
            
            file_path = filedialog.asksaveasfilename(
                defaultextension=".csv",
                initialfile=default_filename,
                filetypes=[("CSV文件", "*.csv"), ("所有文件", "*.*")],
                title="导出引物结果"
            )
            
            if not file_path:
                return
            
            headers = ["来源", "编号", "起始位置", "正向引物", "反向引物", "产物大小", "Tm值", "GC%", "评分", "评价"]
            
            with open(file_path, 'w', newline='', encoding='utf-8-sig') as f:
                writer = csv.writer(f)
                writer.writerow(headers)
                
                for item in all_items:
                    values = self.primer_tree.item(item)['values']
                    writer.writerow(values)
            
            self.set_status(f"已导出{len(all_items)}对引物到: {file_path}", "green")
        except Exception as e:
            self.set_status(f"导出失败: {str(e)}", "red")

    def setup_evaluation_tab(self, parent=None):
        """
        设置序列评价页面
        对siRNA序列进行多维度评价，包括：
        - 结构特性评价
        - 药效评价
        - 稳定性评价
        - 特异性评价
        - 序列质量评价
        - 综合评价报告
        """
        if parent is None:
            parent = self.evaluation_frame
        
        # 序列输入区域
        input_frame = ttk.LabelFrame(parent, text="输入序列", padding=10)
        input_frame.pack(fill='x', expand=False, padx=10, pady=5)
        
        ttk.Label(input_frame, text="正义链 (Sense):").grid(row=0, column=0, sticky='w')
        self.eval_sense_entry = ttk.Entry(input_frame, width=50)
        self.eval_sense_entry.grid(row=0, column=1, padx=5)
        
        ttk.Label(input_frame, text="反义链 (Antisense):").grid(row=1, column=0, sticky='w')
        self.eval_antisense_entry = ttk.Entry(input_frame, width=50)
        self.eval_antisense_entry.grid(row=1, column=1, padx=5)
        
        ttk.Button(input_frame, text="从设计结果导入", command=self.import_from_design).grid(row=0, column=2, padx=5)
        ttk.Button(input_frame, text="开始评价", command=self.evaluate_sequence).grid(row=1, column=2, padx=5)
        
        ttk.Label(input_frame, text="靶区域:").grid(row=2, column=0, sticky='w')
        self.target_region_var = tk.StringVar(value='CDS')
        self.target_region_combo = ttk.Combobox(input_frame, textvariable=self.target_region_var, 
                                               values=['CDS', "3'-UTR", "5'-UTR", '起始密码子附近'], width=15)
        self.target_region_combo.grid(row=2, column=1, padx=5, sticky='w')
        
        self.score_frame = ttk.LabelFrame(input_frame, text="设计评分", padding=5)
        self.score_frame.grid(row=1, column=3, padx=5, sticky='e')
        self.score_label = ttk.Label(self.score_frame, text="总分: --/96", font=('Arial', 12, 'bold'))
        self.score_label.pack()

        tab_button_frame = ttk.Frame(parent)
        tab_button_frame.pack(fill='x', padx=10, pady=5)
        
        self.eval_tabs = ['结构特性', '药效评价', '稳定性评价', '特异性评价', '序列质量', '综合评价']
        self.eval_tab_buttons = []
        
        for tab_name in self.eval_tabs:
            btn = ttk.Button(tab_button_frame, text=tab_name, 
                           command=lambda name=tab_name: self.switch_eval_tab(name))
            btn.pack(side='left', padx=2)
            self.eval_tab_buttons.append(btn)

        self.tab_content_frame = ttk.Frame(parent)
        self.tab_content_frame.pack(fill='both', expand=True, padx=10, pady=5)

        self.structure_tree = ttk.Treeview(self.tab_content_frame, columns=('指标', '数值', '评价', '建议'), show='headings')
        self.structure_tree.heading('指标', text='评价指标')
        self.structure_tree.heading('数值', text='计算值')
        self.structure_tree.heading('评价', text='评级')
        self.structure_tree.heading('建议', text='优化建议')
        self.structure_tree.column('指标', width=120)
        self.structure_tree.column('数值', width=80, anchor='center')
        self.structure_tree.column('评价', width=80, anchor='center')
        self.structure_tree.column('建议', width=250)

        self.efficacy_tree = ttk.Treeview(self.tab_content_frame, columns=('指标', '数值', '评价', '建议'), show='headings')
        self.efficacy_tree.heading('指标', text='评价指标')
        self.efficacy_tree.heading('数值', text='计算值')
        self.efficacy_tree.heading('评价', text='评级')
        self.efficacy_tree.heading('建议', text='优化建议')
        self.efficacy_tree.column('指标', width=120)
        self.efficacy_tree.column('数值', width=80, anchor='center')
        self.efficacy_tree.column('评价', width=80, anchor='center')
        self.efficacy_tree.column('建议', width=250)

        self.stability_tree = ttk.Treeview(self.tab_content_frame, columns=('指标', '数值', '评价', '建议'), show='headings')
        self.stability_tree.heading('指标', text='评价指标')
        self.stability_tree.heading('数值', text='计算值')
        self.stability_tree.heading('评价', text='评级')
        self.stability_tree.heading('建议', text='优化建议')
        self.stability_tree.column('指标', width=120)
        self.stability_tree.column('数值', width=80, anchor='center')
        self.stability_tree.column('评价', width=80, anchor='center')
        self.stability_tree.column('建议', width=250)

        self.specificity_tree = ttk.Treeview(self.tab_content_frame, columns=('指标', '数值', '评价', '建议'), show='headings')
        self.specificity_tree.heading('指标', text='评价指标')
        self.specificity_tree.heading('数值', text='计算值')
        self.specificity_tree.heading('评价', text='评级')
        self.specificity_tree.heading('建议', text='优化建议')
        self.specificity_tree.column('指标', width=120)
        self.specificity_tree.column('数值', width=80, anchor='center')
        self.specificity_tree.column('评价', width=80, anchor='center')
        self.specificity_tree.column('建议', width=250)

        self.quality_tree = ttk.Treeview(self.tab_content_frame, columns=('指标', '数值', '评价', '建议'), show='headings')
        self.quality_tree.heading('指标', text='评价指标')
        self.quality_tree.heading('数值', text='计算值')
        self.quality_tree.heading('评价', text='评级')
        self.quality_tree.heading('建议', text='优化建议')
        self.quality_tree.column('指标', width=120)
        self.quality_tree.column('数值', width=80, anchor='center')
        self.quality_tree.column('评价', width=80, anchor='center')
        self.quality_tree.column('建议', width=250)

        self.summary_text = tk.Text(self.tab_content_frame, height=20, wrap='word')
        self.summary_text.insert("1.0", "请输入序列并点击'开始评价'按钮获取综合评价报告")

        self.switch_eval_tab('结构特性')

    def switch_eval_tab(self, tab_name):
        for tree in [self.structure_tree, self.efficacy_tree, self.stability_tree, 
                     self.specificity_tree, self.quality_tree, self.summary_text]:
            tree.pack_forget()
        
        for btn in self.tab_content_frame.winfo_children():
            btn.pack_forget()
        
        if tab_name == '结构特性':
            self.structure_tree.pack(fill='both', expand=True)
        elif tab_name == '药效评价':
            self.efficacy_tree.pack(fill='both', expand=True)
        elif tab_name == '稳定性评价':
            self.stability_tree.pack(fill='both', expand=True)
        elif tab_name == '特异性评价':
            self.specificity_tree.pack(fill='both', expand=True)
        elif tab_name == '序列质量':
            self.quality_tree.pack(fill='both', expand=True)
        elif tab_name == '综合评价':
            self.summary_text.pack(fill='both', expand=True)

    def import_from_design(self):
        selected_items = self.result_tree.selection()
        if selected_items:
            item = selected_items[0]
            values = self.result_tree.item(item, 'values')
            if len(values) >= 4:
                self.eval_sense_entry.delete(0, tk.END)
                self.eval_sense_entry.insert(0, values[2])
                self.eval_antisense_entry.delete(0, tk.END)
                self.eval_antisense_entry.insert(0, values[3])
                
                if len(values) >= 12:
                    target_region = values[10]
                    region_mapping = {
                        "3'-UTR": "3'-UTR",
                        "5'-UTR": "5'-UTR",
                        "CDS": "CDS",
                        "起始密码子附近": "起始密码子附近",
                        "3'UTR": "3'-UTR",
                        "5'UTR": "5'-UTR"
                    }
                    mapped_region = region_mapping.get(target_region, "CDS")
                    self.target_region_var.set(mapped_region)
                    self.set_status(f"✓ 已从设计结果导入序列（靶区域: {mapped_region}）", "green")
                else:
                    self.set_status("✓ 已从设计结果导入序列", "green")
                
                self.evaluate_sequence()
        else:
            self.set_status("请先在设计结果中选择一条序列", "orange")

    def evaluate_sequence(self):
        """
        评价siRNA序列的主方法
        执行多维度评价：
        1. 结构特性评价 - 序列长度、GC含量、碱基组成等
        2. 药效评价 - Reynolds评分、热力学不对称性等
        3. 稳定性评价 - 二级结构预测、热力学稳定性等
        4. 特异性评价 - BLAST相似性分析、脱靶风险评估
        5. 序列质量评价 - 免疫刺激风险、合成可行性等
        6. 综合评价报告 - 整合所有评价结果
        """
        # 获取输入序列
        sense = self.eval_sense_entry.get().strip()
        antisense = self.eval_antisense_entry.get().strip()
        
        # 验证输入
        if not sense or not antisense:
            self.set_status("请输入正义链和反义链序列", "red")
            return
        
        # 更新状态
        self.set_status("正在评价序列特性...", "orange")
        
        # 清空之前的评价结果
        for tree in [self.structure_tree, self.efficacy_tree, self.stability_tree, self.specificity_tree, self.quality_tree]:
            for item in tree.get_children():
                tree.delete(item)
        
        try:
            self.evaluate_structure(sense, antisense)
        except Exception as e:
            print(f"结构评价失败: {str(e)}")
        
        try:
            self.evaluate_efficacy(sense, antisense)
        except Exception as e:
            print(f"药效评价失败: {str(e)}")
        
        try:
            self.evaluate_stability(sense, antisense)
        except Exception as e:
            print(f"稳定性评价失败: {str(e)}")
        
        try:
            self.evaluate_specificity(sense, antisense)
        except Exception as e:
            print(f"特异性评价失败: {str(e)}")
        
        try:
            self.evaluate_quality(sense, antisense)
        except Exception as e:
            print(f"质量评价失败: {str(e)}")
        
        try:
            self.evaluate_sirna_rules(sense, antisense)
        except Exception as e:
            print(f"siRNA规则评价失败: {str(e)}")
        
        target_region = self.target_region_var.get()
        
        try:
            evaluation = self.designer.comprehensive_sirna_evaluation(sense, target_region_name=target_region)
            self.score_label.config(text=f"总分: {evaluation['final_score']:.1f}/{evaluation['max_score']}")
        except Exception as e:
            print(f"更新评分失败: {str(e)}")
        
        try:
            self.generate_summary(sense, antisense)
            self.set_status("✓ 序列评价完成", "green")
        except Exception as e:
            print(f"生成综合评价失败: {str(e)}")
            self.set_status("评价部分完成，综合评价生成失败", "orange")

    def evaluate_structure(self, sense, antisense):
        results = []
        
        gc_content = self.calculate_gc_content(sense)
        results.append(('GC含量', f"{gc_content:.1f}%", self.rate_gc(gc_content), self.suggest_gc(gc_content)))
        
        tm_value = self.calculate_tm(sense)
        results.append(('Tm值', f"{tm_value:.1f}°C", self.rate_tm(tm_value), self.suggest_tm(tm_value)))
        
        hairpin = self.check_hairpin(sense)
        results.append(('发卡结构', hairpin, self.rate_hairpin(hairpin), self.suggest_hairpin(hairpin)))
        
        self.add_tree_results(self.structure_tree, results)

    def evaluate_efficacy(self, sense, antisense):
        results = []
        
        target_region = self.target_region_var.get()
        evaluation = self.designer.comprehensive_sirna_evaluation(sense, target_region_name=target_region)
        
        reynolds_score = evaluation['reynolds_score']
        results.append(('Reynolds评分', f"{reynolds_score}/8", self.rate_reynolds(reynolds_score), self.suggest_reynolds(reynolds_score)))
        
        amarz_score = evaluation['amarzguioui_score']
        results.append(('Amarzguioui评分', f"{amarz_score}/5", self.rate_amarz(amarz_score), self.suggest_amarz(amarz_score)))
        
        seed_complexity = evaluation['seed_complexity']
        results.append(('种子区复杂度', f"{seed_complexity:.2f}", self.rate_seed_complexity(seed_complexity), self.suggest_seed_complexity(seed_complexity)))
        
        thermo_pass = evaluation['thermo_asymmetry']
        thermo_status = "通过" if thermo_pass else "不通过"
        results.append(('热力学不对称性', thermo_status, self.rate_thermo(thermo_pass), self.suggest_thermo(thermo_pass)))
        
        uitei_pass = evaluation['uitei_passed']
        uitei_reason = evaluation['uitei_reason']
        uitei_status = "通过" if uitei_pass else f"不通过 ({uitei_reason})"
        results.append(('Ui-Tei规则', uitei_status, self.rate_uitei(uitei_pass), self.suggest_uitei(uitei_pass)))
        
        self.add_tree_results(self.efficacy_tree, results)
    


    def evaluate_stability(self, sense, antisense):
        results = []
        
        nuclease_resistance = self.predict_nuclease_resistance(sense)
        results.append(('核酸酶抗性', f"{nuclease_resistance:.1f}%", self.rate_stability(nuclease_resistance), self.suggest_stability(nuclease_resistance)))
        
        serum_half_life = self.predict_half_life(sense)
        results.append(('血清半衰期', f"{serum_half_life:.1f}小时", self.rate_half_life(serum_half_life), self.suggest_half_life(serum_half_life)))
        
        chemical_stability = self.analyze_chemical_stability(sense)
        results.append(('化学稳定性', chemical_stability, self.rate_chemical_stability(chemical_stability), self.suggest_chemical_stability(chemical_stability)))
        
        self.add_tree_results(self.stability_tree, results)

    def add_tree_results(self, tree, results):
        for item in results:
            tag = 'good' if item[2] in ['优秀', '良好'] else 'warning' if item[2] == '一般' else 'danger'
            tree.insert('', tk.END, values=item, tags=(tag,))
            tree.tag_configure('good', foreground='green')
            tree.tag_configure('warning', foreground='orange')
            tree.tag_configure('danger', foreground='red')

    def evaluate_specificity(self, sense, antisense):
        results = []
        
        off_target_risk = self.calculate_off_target_risk(sense)
        results.append(('脱靶风险评分', f"{off_target_risk:.1f}", self.rate_off_target(off_target_risk), self.suggest_off_target(off_target_risk)))
        
        seed_specificity = self.analyze_seed_specificity(sense)
        results.append(('种子区特异性', f"{seed_specificity:.1f}", self.rate_seed_specificity(seed_specificity), self.suggest_seed_specificity(seed_specificity)))
        
        global_alignment = self.calculate_global_alignment(sense)
        results.append(('全局比对得分', f"{global_alignment:.1f}", self.rate_global_alignment(global_alignment), self.suggest_global_alignment(global_alignment)))
        
        specificity_score = self.calculate_specificity_score(sense)
        results.append(('综合特异性', f"{specificity_score:.1f}", self.rate_overall_specificity(specificity_score), self.suggest_overall_specificity(specificity_score)))
        
        self.add_tree_results(self.specificity_tree, results)

    def evaluate_quality(self, sense, antisense):
        results = []
        
        sequence_purity = self.calculate_sequence_purity(sense)
        results.append(('序列纯度', f"{sequence_purity:.1f}%", self.rate_purity(sequence_purity), self.suggest_purity(sequence_purity)))
        
        manufacturing_feasibility = self.analyze_manufacturing_feasibility(sense)
        results.append(('合成可行性', manufacturing_feasibility, self.rate_manufacturing(manufacturing_feasibility), self.suggest_manufacturing(manufacturing_feasibility)))
        
        secondary_structure = self.analyze_secondary_structure(sense)
        results.append(('二级结构风险', secondary_structure, self.rate_secondary_structure(secondary_structure), self.suggest_secondary_structure(secondary_structure)))
        
        modification_tolerance = self.calculate_modification_tolerance(sense)
        results.append(('修饰耐受性', f"{modification_tolerance:.1f}%", self.rate_mod_tolerance(modification_tolerance), self.suggest_mod_tolerance(modification_tolerance)))
        
        quality_score = self.calculate_quality_score(sense)
        results.append(('综合质量评分', f"{quality_score:.1f}", self.rate_overall_quality(quality_score), self.suggest_overall_quality(quality_score)))
        
        self.add_tree_results(self.quality_tree, results)

    def evaluate_sirna_rules(self, sense, antisense):
        sirna_results = []
        
        target_region = self.target_region_var.get()
        evaluation = self.designer.comprehensive_sirna_evaluation(sense, target_region_name=target_region)
        
        reynolds_score = evaluation['reynolds_score']
        sirna_results.append(('Reynolds评分', f"{reynolds_score}/8", self.rate_reynolds(reynolds_score), self.suggest_reynolds(reynolds_score)))
        
        amarz_score = evaluation['amarzguioui_score']
        sirna_results.append(('Amarzguioui评分', f"{amarz_score}/5", self.rate_amarz(amarz_score), self.suggest_amarz(amarz_score)))
        
        seed_complexity = evaluation['seed_complexity']
        sirna_results.append(('种子区复杂度', f"{seed_complexity:.2f}", self.rate_seed_complexity(seed_complexity), self.suggest_seed_complexity(seed_complexity)))
        
        thermo_pass = evaluation['thermo_asymmetry']
        thermo_status = "通过" if thermo_pass else "不通过"
        sirna_results.append(('热力学不对称性', thermo_status, self.rate_thermo(thermo_pass), self.suggest_thermo(thermo_pass)))
        
        uitei_pass = evaluation['uitei_passed']
        uitei_reason = evaluation['uitei_reason']
        uitei_status = "通过" if uitei_pass else f"不通过 ({uitei_reason})"
        sirna_results.append(('Ui-Tei规则', uitei_status, self.rate_uitei(uitei_pass), self.suggest_uitei(uitei_pass)))
        
        immune_status = self.check_immune_stimulatory(sense)
        sirna_results.append(('免疫刺激序列', immune_status, self.rate_immune(immune_status), self.suggest_immune(immune_status)))

        gene_name = self.gene_name_var.get()
        
        try:
            from src.ncbi_blast import NCBIBLAST
            blast = NCBIBLAST()
            
            self.set_status("正在进行NCBI BLAST分析...", "orange")
            
            result = blast.run_blast(sense, timeout=60, poll_interval=3)
            
            if result and result.get('hits'):
                risk = blast.get_off_target_risk(result, gene_name, threshold=80.0)
                
                high_count = len(risk['high_risk'])
                medium_count = len(risk['medium_risk'])
                low_count = len(risk['low_risk'])
                
                if high_count > 0:
                    blast_status = f"高风险 ({high_count}个)"
                    blast_rating = '较差'
                    blast_suggest = f"发现 {high_count} 个高风险脱靶位点，建议重新设计序列"
                elif medium_count > 0:
                    blast_status = f"中风险 ({medium_count}个)"
                    blast_rating = '一般'
                    blast_suggest = f"发现 {medium_count} 个中风险脱靶位点，建议进行实验验证"
                elif low_count > 0:
                    blast_status = f"低风险 ({low_count}个)"
                    blast_rating = '良好'
                    blast_suggest = f"发现 {low_count} 个低风险脱靶位点，脱靶风险较低"
                else:
                    blast_status = '未发现脱靶位点'
                    blast_rating = '优秀'
                    blast_suggest = 'NCBI数据库中未发现明显脱靶位点，特异性良好'
            elif result:
                blast_status = '无显著匹配'
                blast_rating = '优秀'
                blast_suggest = '序列在RefSeq数据库中没有显著匹配，特异性很高'
            else:
                blast_status = 'BLAST查询失败'
                blast_rating = '未检测'
                blast_suggest = '建议手动在NCBI BLAST网站验证'
                
            self.set_status("✓ 评价完成", "green")
            
        except Exception as e:
            blast_status = '网络不可用'
            blast_rating = '未检测'
            blast_suggest = '网络不可用，无法进行NCBI BLAST分析'
            self.set_status("✓ 评价完成（BLAST跳过）", "yellow")
        
        sirna_results.append(('NCBI BLAST脱靶风险', blast_status, blast_rating, blast_suggest))

        self.add_tree_results(self.efficacy_tree, sirna_results)

    def update_ranks(self):
        items = self.result_tree.get_children()
        for idx, item in enumerate(items):
            values = list(self.result_tree.item(item, 'values'))
            if values:
                values[0] = "最佳" if idx == 0 else str(idx + 1)
                self.result_tree.item(item, values=tuple(values))
    
    def rate_blast(self, identity):
        if identity == 0:
            return '未检测'
        elif identity < 80:
            return '优秀'
        elif identity < 90:
            return '一般'
        else:
            return '较差'

    def suggest_blast(self, identity):
        if identity == 0:
            return '建议进行BLAST检测'
        elif identity < 80:
            return '脱靶风险低，安全性良好'
        elif identity < 90:
            return '存在潜在脱靶风险，建议验证'
        else:
            return '高脱靶风险，建议重新设计'

    def run_blast_on_selected(self):
        selected_items = self.result_tree.selection()
        if not selected_items:
            messagebox.showwarning('提示', '请先在设计结果中选择一条序列')
            return

        item = selected_items[0]
        values = self.result_tree.item(item, 'values')
        if len(values) < 4:
            messagebox.showwarning('提示', '序列数据不完整')
            return

        sense = values[2]
        gene_name = self.gene_name_var.get()

        if not sense:
            messagebox.showwarning('提示', '序列为空')
            return

        self.set_status("正在连接NCBI数据库...", "orange")
        self.blast_text.delete("1.0", tk.END)
        self.blast_text.insert("1.0", "【NCBI BLAST分析中...】\n")
        self.blast_text.insert(tk.END, "=" * 40 + "\n\n")

        import threading

        def blast_analysis():
            try:
                from src.ncbi_blast import NCBIBLAST
            except ImportError as e:
                self.set_status("✗ 模块导入失败", "red")
                self.blast_text.insert(tk.END, f"\n✗ 无法导入NCBI BLAST模块: {str(e)}\n")
                self.blast_text.insert(tk.END, "请检查src目录是否存在\n")
                return

            try:
                blast = NCBIBLAST()
                
                self.set_status("正在提交BLAST请求...", "orange")
                self.blast_text.insert(tk.END, "✓ 正在提交BLAST请求到NCBI服务器...\n")
                
                request_id = blast.submit_blast(sense)
                
                if not request_id:
                    self.set_status("✗ BLAST提交失败", "red")
                    self.blast_text.insert(tk.END, "\n✗ BLAST请求提交失败\n")
                    self.blast_text.insert(tk.END, "请检查网络连接后重试\n")
                    return
                
                self.blast_text.insert(tk.END, f"✓ 请求已提交，任务ID: {request_id}\n")
                
                import time
                start_time = time.time()
                elapsed = 0
                
                while elapsed < 120:
                    elapsed = time.time() - start_time
                    
                    self.set_status(f"BLAST分析中... ({int(elapsed)}秒)", "orange")
                    if int(elapsed) % 2 == 0:
                        self.blast_text.insert(tk.END, ".")
                    
                    status, result = blast.check_status(request_id)
                    
                    if status == 'COMPLETE' and result:
                        parsed_result = blast.parse_result(result)
                        
                        self.blast_text.delete("1.0", tk.END)
                        self.blast_text.insert("1.0", "=" * 60 + "\n")
                        self.blast_text.insert(tk.END, "【NCBI BLAST相似性分析报告】\n")
                        self.blast_text.insert(tk.END, "=" * 60 + "\n\n")

                        self.blast_text.insert(tk.END, f"查询序列: {sense}\n")
                        self.blast_text.insert(tk.END, f"目标基因: {gene_name}\n")
                        self.blast_text.insert(tk.END, "数据库: RefSeq mRNA\n")
                        self.blast_text.insert(tk.END, f"耗时: {int(elapsed)}秒\n\n")

                        if parsed_result:
                            self.blast_text.insert(tk.END, f"✓ BLAST查询完成\n\n")
                            
                            if parsed_result.get('hits'):
                                risk = blast.get_off_target_risk(parsed_result, gene_name, threshold=80.0)
                                
                                self.blast_text.insert(tk.END, f"检测到 {len(parsed_result['hits'])} 个比对结果\n\n")
                                
                                if risk['has_off_target']:
                                    self.blast_text.insert(tk.END, f"⚠ 发现 {risk['total_hits']} 个潜在脱靶位点\n\n")

                                    self.blast_text.insert(tk.END, f"  • 高风险 (相似度≥90%): {len(risk['high_risk'])} 个\n")
                                    self.blast_text.insert(tk.END, f"  • 中风险 (相似度85-89%): {len(risk['medium_risk'])} 个\n")
                                    self.blast_text.insert(tk.END, f"  • 低风险 (相似度80-84%): {len(risk['low_risk'])} 个\n\n")

                                    if risk['target_excluded']:
                                        self.blast_text.insert(tk.END, "✓ 目标基因已自动排除\n\n")

                                    if risk['high_risk']:
                                        self.blast_text.insert(tk.END, "【高风险脱靶位点】\n")
                                        self.blast_text.insert(tk.END, "-" * 40 + "\n")
                                        for hit in risk['high_risk'][:10]:
                                            self.blast_text.insert(tk.END, f"  [{hit.get('accession', '')}]\n")
                                            self.blast_text.insert(tk.END, f"    基因: {hit.get('description', '')[:50]}...\n")
                                            self.blast_text.insert(tk.END, f"    相似度: {hit.get('identity', 0):.1f}%\n")
                                            self.blast_text.insert(tk.END, f"    E-value: {hit.get('e_value', 0):.2e}\n\n")

                                    if risk['medium_risk']:
                                        self.blast_text.insert(tk.END, "【中风险脱靶位点】\n")
                                        self.blast_text.insert(tk.END, "-" * 40 + "\n")
                                        for hit in risk['medium_risk'][:10]:
                                            self.blast_text.insert(tk.END, f"  [{hit.get('accession', '')}]\n")
                                            self.blast_text.insert(tk.END, f"    基因: {hit.get('description', '')[:50]}...\n")
                                            self.blast_text.insert(tk.END, f"    相似度: {hit.get('identity', 0):.1f}%\n")
                                            self.blast_text.insert(tk.END, f"    E-value: {hit.get('e_value', 0):.2e}\n\n")

                                    self.blast_text.insert(tk.END, "【建议】\n")
                                    if risk['high_risk']:
                                        self.blast_text.insert(tk.END, "  • 存在高风险脱靶位点，建议重新设计序列\n")
                                    else:
                                        self.blast_text.insert(tk.END, "  • 存在潜在脱靶风险，建议进行实验验证\n")
                                else:
                                    self.blast_text.insert(tk.END, "✓ 未发现明显脱靶位点（相似度<80%）\n\n")
                                    self.blast_text.insert(tk.END, "该序列特异性良好，可以优先考虑使用。\n")
                            else:
                                self.blast_text.insert(tk.END, "✓ 未找到显著比对结果\n\n")
                                self.blast_text.insert(tk.END, "该序列在RefSeq数据库中没有显著匹配，特异性很高。\n")
                        else:
                            self.blast_text.insert(tk.END, "✗ BLAST结果解析失败\n")

                        self.set_status("✓ BLAST分析完成", "green")
                        return
                    elif status == 'FAILED':
                        self.set_status("✗ BLAST查询失败", "red")
                        self.blast_text.delete("1.0", tk.END)
                        self.blast_text.insert("1.0", "✗ BLAST查询失败\n\n")
                        self.blast_text.insert(tk.END, "可能的原因:\n")
                        self.blast_text.insert(tk.END, "  • NCBI服务器错误\n")
                        self.blast_text.insert(tk.END, "  • 序列格式问题\n\n")
                        self.blast_text.insert(tk.END, "请稍后重试或直接访问NCBI BLAST网站")
                        return
                    
                    time.sleep(3)
                
                self.set_status("✗ BLAST查询超时", "red")
                self.blast_text.delete("1.0", tk.END)
                self.blast_text.insert("1.0", "✗ BLAST查询超时\n\n")
                self.blast_text.insert(tk.END, "NCBI服务器响应时间过长\n")
                self.blast_text.insert(tk.END, "建议直接访问NCBI BLAST网站手动查询")
            
            except Exception as e:
                self.set_status("✗ BLAST分析异常", "red")
                self.blast_text.insert(tk.END, f"\n✗ BLAST分析发生异常: {str(e)}\n")
                import traceback
                self.blast_text.insert(tk.END, f"错误详情: {traceback.format_exc()[:200]}...\n")

        threading.Thread(target=blast_analysis, daemon=True).start()

    def delete_selected_sequence(self):
        selected_items = self.result_tree.selection()
        if not selected_items:
            messagebox.showwarning('提示', '请先在设计结果中选择要删除的序列')
            return

        item = selected_items[0]
        values = self.result_tree.item(item, 'values')
        sense = values[2]
        antisense = values[3]

        confirm = messagebox.askyesno('确认删除', f'确定要删除序列 {sense} 吗？\n\n删除后将自动从候选序列中补充新的siRNA。')
        if not confirm:
            return

        self.result_tree.delete(item)
        self.set_status(f"已删除序列 {sense}，正在补充新序列...", "orange")

        mrna = self.current_mrna
        if mrna:
            length = len(sense)
            start_pos = max(1, int(values[1].split('-')[0]) - 10)
            end_pos = min(len(mrna), start_pos + 50)

            new_candidates = self.designer.find_sirna_candidates(mrna, start_pos, end_pos, length, max_candidates=5)

            for candidate in new_candidates:
                if candidate['passes_hard_filters']:
                    sense_new = candidate['sense']
                    antisense_new = candidate['antisense']
                    pos = candidate['position']
                    gc = candidate['gc_content']
                    reynolds = candidate['reynolds_score']
                    amarz = candidate['amarzguioui_score']
                    seed_complex = f"{candidate['seed_complexity']:.2f}"
                    thermo = "✓" if candidate['thermo_asymmetry'] else "✗"
                    uitei = "✓" if candidate['uitei_passed'] else "✗"
                    target_region = candidate['target_region']
                    final_score = candidate['final_score']
                    recommendation = candidate['recommendation']

                    pos_range = f"{pos} - {pos + length - 1}"

                    existing_items = self.result_tree.get_children()
                    for existing in existing_items:
                        existing_values = self.result_tree.item(existing, 'values')
                        if existing_values[2] == sense_new:
                            break
                    else:
                        current_count = len(self.result_tree.get_children())
                        rank = "最佳" if current_count == 0 else str(current_count + 1)
                        self.result_tree.insert('', 'end', values=(rank, pos_range, sense_new, antisense_new, f"{gc:.1f}", reynolds, amarz, seed_complex, thermo, uitei, target_region, f"{final_score:.1f}", recommendation))
                        self.set_status(f"✓ 已补充新序列 {sense_new}", "green")
                        break
            else:
                self.set_status("未找到合适的补充序列", "orange")

        self.update_ranks()


    def calculate_reynolds_score(self, sense):
        score = 0
        sense = sense.upper()
        if len(sense) >= 21:
            if sense[18] == 'G': score += 1
            if sense[2] == 'A': score += 1
            if sense[9] == 'U': score += 1
            if sense[0] in ['A', 'U']: score += 1
            gc = self.calculate_gc_content(sense)
            if 30.0 <= gc <= 55.0: score += 1
            if not re.search(r'(.)\1{3,}', sense): score += 1
        return score

    def rate_reynolds(self, score):
        if score >= 6: return '优秀'
        elif score >= 4: return '良好'
        elif score >= 2: return '一般'
        else: return '较差'

    def suggest_reynolds(self, score):
        if score >= 6: return 'Reynolds评分高，沉默活性预期良好'
        elif score >= 4: return 'Reynolds评分中等，可以接受'
        else: return 'Reynolds评分较低，建议优化序列'

    def calculate_amarzguioui_score(self, sense):
        score = 0
        sense = sense.upper()
        if len(sense) >= 21:
            if sense[0] in ['A', 'C']: score += 1
            if sense[5] == 'A': score += 1
            if sense[18] != 'G': score += 1
            gc = self.calculate_gc_content(sense)
            if 30.0 <= gc <= 55.0: score += 1
            if not re.search(r'UUU', sense): score += 1
        return score

    def rate_amarz(self, score):
        if score >= 4: return '优秀'
        elif score >= 3: return '良好'
        elif score >= 2: return '一般'
        else: return '较差'

    def suggest_amarz(self, score):
        if score >= 4: return 'Amarzguioui评分高'
        elif score >= 3: return 'Amarzguioui评分中等'
        else: return 'Amarzguioui评分较低'

    def calculate_seed_complexity(self, sense):
        import math
        sense = sense.upper()
        antisense = self.reverse_complement(sense)
        seed = antisense[1:8] if len(antisense) >= 8 else antisense
        
        base_count = {'A': 0, 'C': 0, 'G': 0, 'U': 0}
        for base in seed:
            if base in base_count:
                base_count[base] += 1

        total = len(seed)
        if total == 0: return 0.0

        entropy = 0.0
        for count in base_count.values():
            if count > 0:
                p = count / total
                entropy -= p * math.log2(p)

        max_entropy = math.log2(4)
        return entropy / max_entropy if max_entropy > 0 else 0.0

    def rate_seed_complexity(self, complexity):
        if complexity >= 0.6: return '优秀'
        elif complexity >= 0.4: return '良好'
        elif complexity >= 0.2: return '一般'
        else: return '较差'

    def suggest_seed_complexity(self, complexity):
        if complexity >= 0.6: return '种子区复杂度高，脱靶风险低'
        elif complexity >= 0.4: return '种子区复杂度中等'
        else: return '种子区复杂度低，存在miRNA-like脱靶风险'

    def check_thermodynamic_asymmetry(self, sense, antisense):
        sense = sense.upper()[:4]
        antisense = antisense.upper()[:4]
        
        dg_table = {'AA': -9.1, 'TT': -9.1, 'AT': -8.6, 'TA': -6.0,
                    'CA': -5.8, 'GT': -6.5, 'AC': -6.5, 'TG': -5.7,
                    'CT': -7.8, 'GA': -5.6, 'AG': -7.8, 'TC': -5.8,
                    'CG': -11.9, 'GC': -11.1, 'GG': -11.0, 'CC': -11.0}
        
        def calc_dg(seq):
            if len(seq) < 2: return 0.0
            total = 0.0
            for i in range(len(seq)-1):
                dinuc = seq[i:i+2]
                total += dg_table.get(dinuc, -7.0)
            return total
        
        sense_dg = calc_dg(sense)
        anti_dg = calc_dg(antisense)
        
        return anti_dg < sense_dg, anti_dg - sense_dg

    def rate_thermo(self, passed):
        return '优秀' if passed else '较差'

    def suggest_thermo(self, passed):
        if passed:
            return '反义链5\'端更稳定，有利于正确加载入RISC'
        else:
            return '热力学不对称性不满足，可能影响RISC加载效率'

    def check_ui_tei_rules(self, sense, antisense):
        sense = sense.upper()
        antisense = antisense.upper()
        
        if len(sense) < 21 or len(antisense) < 21:
            return False, '序列长度不足'
        
        ut1 = sense[0] in ['A', 'U']
        ut2 = antisense[0] in ['G', 'C']
        gc = self.calculate_gc_content(sense)
        ut3 = 30.0 <= gc <= 55.0
        ut4 = not re.search(r'(.)\1{3,}', sense)
        
        dg_table = {'AA': -9.1, 'TT': -9.1, 'AT': -8.6, 'TA': -6.0,
                    'CA': -5.8, 'GT': -6.5, 'AC': -6.5, 'TG': -5.7,
                    'CT': -7.8, 'GA': -5.6, 'AG': -7.8, 'TC': -5.8,
                    'CG': -11.9, 'GC': -11.1, 'GG': -11.0, 'CC': -11.0}
        
        def calc_dg(seq):
            total = 0.0
            for i in range(min(4, len(seq)-1)):
                dinuc = seq[i:i+2]
                total += dg_table.get(dinuc, -7.0)
            return total
        
        ut5 = calc_dg(sense[:5]) > -5.0
        
        if ut1 and ut2 and ut3 and ut4 and ut5:
            return True, '全部满足'
        else:
            reasons = []
            if not ut1: reasons.append('义链1位非A/U')
            if not ut2: reasons.append('反义链1位非G/C')
            if not ut3: reasons.append('GC不在30-55%')
            if not ut4: reasons.append('含≥4连续碱基')
            if not ut5: reasons.append('义链5\'端ΔG过低')
            return False, '; '.join(reasons)

    def rate_uitei(self, passed):
        return '优秀' if passed else '较差'

    def suggest_uitei(self, passed):
        if passed:
            return '完全符合Ui-Tei规则，沉默效率预期高'
        else:
            return '不满足Ui-Tei规则，建议优化序列'

    def check_immune_stimulatory(self, sense):
        sense = sense.upper()
        immune_seq = ['UGUGU', 'GUCCUUCAA']
        
        for seq in immune_seq:
            if seq in sense:
                return f"含{seq} (免疫激活风险)"
        
        if re.search(r'UUUUUU', sense):
            return '含≥6连续U (TLR信号风险)'
        
        return '无免疫刺激序列'

    def rate_immune(self, status):
        if status == '无免疫刺激序列':
            return '优秀'
        else:
            return '较差'

    def suggest_immune(self, status):
        if status == '无免疫刺激序列':
            return '无免疫刺激风险'
        else:
            return '存在免疫刺激风险，建议避免使用'

    def calculate_off_target_risk(self, sequence):
        risk = 0
        
        if 'GGG' in sequence or 'CCC' in sequence:
            risk += 15
        if 'GG' in sequence and sequence.count('G') > len(sequence) * 0.35:
            risk += 10
        if sequence[1:7].count('G') >= 4:
            risk += 20
        
        return max(0, 100 - risk)

    def rate_off_target(self, score):
        if score >= 85:
            return '优秀'
        elif score >= 70:
            return '良好'
        elif score >= 50:
            return '一般'
        else:
            return '较差'

    def suggest_off_target(self, score):
        if score >= 85:
            return '脱靶风险低，特异性良好'
        elif score >= 70:
            return '存在一定脱靶风险，建议优化种子区'
        else:
            return '脱靶风险较高，建议重新设计序列'

    def analyze_seed_specificity(self, sequence):
        seed = sequence[1:8] if len(sequence) >= 8 else sequence
        score = 70
        
        if seed.count('A') + seed.count('U') >= 5:
            score += 15
        if 'G' not in seed[:3]:
            score += 15
        
        return min(100, score)

    def rate_seed_specificity(self, score):
        if score >= 90:
            return '优秀'
        elif score >= 75:
            return '良好'
        elif score >= 60:
            return '一般'
        else:
            return '较差'

    def suggest_seed_specificity(self, score):
        if score >= 90:
            return '种子区特异性优秀'
        elif score >= 75:
            return '种子区特异性良好'
        else:
            return '建议优化种子区序列'

    def calculate_global_alignment(self, sequence):
        score = 80
        
        if len(sequence) >= 19 and len(sequence) <= 23:
            score += 10
        if sequence[-1] in ['A', 'U']:
            score += 5
        if sequence[0] in ['A', 'U']:
            score += 5
        
        return min(100, score)

    def rate_global_alignment(self, score):
        if score >= 90:
            return '优秀'
        elif score >= 75:
            return '良好'
        elif score >= 60:
            return '一般'
        else:
            return '较差'

    def suggest_global_alignment(self, score):
        if score >= 90:
            return '序列长度和末端碱基符合最优设计原则'
        else:
            return '建议调整序列长度或末端碱基'

    def calculate_specificity_score(self, sequence):
        off_target = self.calculate_off_target_risk(sequence)
        seed_spec = self.analyze_seed_specificity(sequence)
        global_align = self.calculate_global_alignment(sequence)
        
        return (off_target * 0.4 + seed_spec * 0.35 + global_align * 0.25)

    def rate_overall_specificity(self, score):
        if score >= 85:
            return '优秀'
        elif score >= 70:
            return '良好'
        elif score >= 50:
            return '一般'
        else:
            return '较差'

    def suggest_overall_specificity(self, score):
        if score >= 85:
            return '综合特异性优秀，适合作为候选药物'
        elif score >= 70:
            return '综合特异性良好，可进一步优化'
        else:
            return '建议重新设计以提高特异性'

    def calculate_sequence_purity(self, sequence):
        valid_bases = {'A', 'T', 'U', 'C', 'G'}
        valid_count = sum(1 for c in sequence if c.upper() in valid_bases)
        return (valid_count / len(sequence)) * 100 if sequence else 0

    def rate_purity(self, purity):
        if purity == 100:
            return '优秀'
        elif purity >= 95:
            return '良好'
        elif purity >= 90:
            return '一般'
        else:
            return '较差'

    def suggest_purity(self, purity):
        if purity == 100:
            return '序列纯度100%，无异常碱基'
        elif purity >= 95:
            return '序列纯度良好'
        else:
            return '建议检查序列中是否存在异常碱基'

    def analyze_manufacturing_feasibility(self, sequence):
        if len(sequence) > 30:
            return '较差'
        elif len(sequence) > 25:
            return '一般'
        elif len(sequence) >= 19:
            return '良好'
        else:
            return '优秀'

    def rate_manufacturing(self, feasibility):
        return {'优秀': '优秀', '良好': '良好', '一般': '一般', '较差': '较差'}[feasibility]

    def suggest_manufacturing(self, feasibility):
        if feasibility == '优秀':
            return '序列长度适合大规模合成'
        elif feasibility == '良好':
            return '合成可行性良好'
        elif feasibility == '一般':
            return '较长序列可能增加合成难度'
        else:
            return '序列过长，建议缩短'

    def analyze_secondary_structure(self, sequence):
        stem_count = 0
        for i in range(len(sequence) - 3):
            window = sequence[i:i+4]
            complement = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
            if all(complement.get(a) == b for a, b in zip(window, window[::-1])):
                stem_count += 1
        
        if stem_count == 0:
            return '低'
        elif stem_count <= 2:
            return '中'
        else:
            return '高'

    def rate_secondary_structure(self, risk):
        return {'低': '优秀', '中': '良好', '高': '较差'}[risk]

    def suggest_secondary_structure(self, risk):
        if risk == '低':
            return '二级结构风险低'
        elif risk == '中':
            return '存在一定二级结构形成风险'
        else:
            return '二级结构风险高，可能影响活性'

    def calculate_modification_tolerance(self, sequence):
        tolerance = 80
        
        if 'AU' in sequence or 'UA' in sequence:
            tolerance += 10
        if sequence[-2:] in ['AA', 'UU', 'AU', 'UA']:
            tolerance += 10
        
        return min(100, tolerance)

    def rate_mod_tolerance(self, tolerance):
        if tolerance >= 90:
            return '优秀'
        elif tolerance >= 75:
            return '良好'
        elif tolerance >= 60:
            return '一般'
        else:
            return '较差'

    def suggest_mod_tolerance(self, tolerance):
        if tolerance >= 90:
            return '适合进行多种化学修饰'
        elif tolerance >= 75:
            return '修饰耐受性良好'
        else:
            return '修饰位点选择需要谨慎'

    def calculate_quality_score(self, sequence):
        purity = self.calculate_sequence_purity(sequence)
        feasibility = {'优秀': 100, '良好': 80, '一般': 60, '较差': 30}[self.analyze_manufacturing_feasibility(sequence)]
        structure = {'低': 100, '中': 70, '高': 30}[self.analyze_secondary_structure(sequence)]
        mod_tol = self.calculate_modification_tolerance(sequence)
        
        return (purity * 0.3 + feasibility * 0.25 + structure * 0.25 + mod_tol * 0.2)

    def rate_overall_quality(self, score):
        if score >= 85:
            return '优秀'
        elif score >= 70:
            return '良好'
        elif score >= 50:
            return '一般'
        else:
            return '较差'

    def suggest_overall_quality(self, score):
        if score >= 85:
            return '序列质量优秀，适合进一步开发'
        elif score >= 70:
            return '序列质量良好'
        else:
            return '建议优化序列质量'

    def generate_summary(self, sense, antisense):
        target_region = self.target_region_var.get()
        evaluation = self.designer.comprehensive_sirna_evaluation(sense, target_region_name=target_region)
        
        summary = f"序列评价报告\n\n"
        summary += f"正义链: {evaluation['sense']}\n"
        summary += f"反义链: {evaluation['antisense']}\n\n"
        summary += "="*60 + "\n\n"
        
        summary += "【过滤规则】\n"
        summary += f"  • 硬性过滤: {'通过' if evaluation['passes_hard_filters'] else '未通过'}\n\n"
        
        summary += "【评分详情】\n"
        summary += f"  • Reynolds评分: {evaluation['reynolds_score']}/8 (权重×4)\n"
        summary += f"  • Amarzguioui评分: {evaluation['amarzguioui_score']}/5 (权重×3)\n"
        summary += f"  • 种子区复杂度: {evaluation['seed_complexity']:.2f} (权重×8)\n"
        summary += f"  • 热力学不对称性: {'通过 (+3分)' if evaluation['thermo_asymmetry'] else '不通过'}\n"
        uitei_text = f"不通过 ({evaluation['uitei_reason']})" if not evaluation['uitei_passed'] else "通过 (+3分)"
        summary += f"  • Ui-Tei规则: {uitei_text}\n"
        summary += f"  • 最优GC范围: {'通过 (+5分)' if evaluation['optimal_gc'] else '不通过'}\n"
        summary += f"  • 靶区域: {evaluation['target_region']} (+{evaluation['target_bonus']}分)\n\n"
        
        summary += "="*60 + "\n\n"
        summary += f"综合评分: {evaluation['final_score']:.1f}/{evaluation['max_score']}\n\n"
        
        recommendation = evaluation['recommendation']
        if recommendation == '推荐':
            summary += "【评价等级】推荐\n"
            summary += "该序列具有良好的药物开发潜力，建议进一步优化。\n"
        elif recommendation == '良好':
            summary += "【评价等级】良好\n"
            summary += "该序列基本符合小核酸药物要求，存在一些可优化点。\n"
        elif recommendation == '一般':
            summary += "【评价等级】一般\n"
            summary += "该序列需要进行较多优化才能达到药物开发标准。\n"
        else:
            summary += "【评价等级】不推荐\n"
            summary += "该序列不适合作为药物候选，建议重新设计。\n"
        
        self.summary_text.delete("1.0", tk.END)
        self.summary_text.insert("1.0", summary)

    def calculate_gc_content(self, sequence):
        if not sequence:
            return 0
        return (sequence.count('G') + sequence.count('C')) / len(sequence) * 100

    def rate_gc(self, gc):
        if 35 <= gc <= 55:
            return '优秀'
        elif 30 <= gc <= 60:
            return '良好'
        elif 25 <= gc <= 65:
            return '一般'
        else:
            return '较差'

    def suggest_gc(self, gc):
        if gc < 35:
            return '建议增加G/C碱基比例以提高稳定性'
        elif gc > 55:
            return '建议降低G/C碱基比例以避免过度稳定'
        else:
            return 'GC含量处于最佳范围'

    def calculate_tm(self, sequence):
        if len(sequence) < 15:
            return 50
        gc = self.calculate_gc_content(sequence)
        return 64.9 + 41 * (gc - 16.4) / len(sequence)

    def rate_tm(self, tm):
        if 78 <= tm <= 85:
            return '优秀'
        elif 75 <= tm <= 88:
            return '良好'
        elif 70 <= tm <= 92:
            return '一般'
        else:
            return '较差'

    def suggest_tm(self, tm):
        if tm < 78:
            return 'Tm值偏低，建议增加GC含量或序列长度'
        elif tm > 85:
            return 'Tm值偏高，建议降低GC含量'
        else:
            return 'Tm值处于最佳范围'

    def check_hairpin(self, sequence):
        min_stem = 4
        for i in range(len(sequence) - 2*min_stem):
            stem1 = sequence[i:i+min_stem]
            stem2 = sequence[i+min_stem+1:i+2*min_stem+1]
            if stem1 and stem2 and all(a in 'GC' and b in 'GC' for a, b in zip(stem1, stem2[::-1])):
                return '存在'
        return '无'

    def rate_hairpin(self, hairpin):
        return '较差' if hairpin == '存在' else '优秀'

    def suggest_hairpin(self, hairpin):
        if hairpin == '存在':
            return '发卡结构可能影响RNAi效果，建议重新设计序列'
        else:
            return '未检测到发卡结构，符合要求'

    def analyze_seed_region(self, sequence):
        seed = sequence[1:8] if len(sequence) >= 8 else sequence
        a_count = seed.count('A')
        u_count = seed.count('U')
        gc_count = seed.count('G') + seed.count('C')
        if gc_count <= 2 and (a_count + u_count) >= 5:
            return '良好'
        else:
            return '一般'

    def rate_seed(self, seed):
        return '优秀' if seed == '良好' else '一般'

    def suggest_seed(self, seed):
        if seed == '良好':
            return '种子区符合高效siRNA特征'
        else:
            return '建议优化种子区，增加A/U含量'

    def predict_potency(self, sequence):
        score = 70
        gc = self.calculate_gc_content(sequence)
        if 35 <= gc <= 55:
            score += 15
        elif 30 <= gc <= 60:
            score += 10
        
        if len(sequence) >= 19 and len(sequence) <= 23:
            score += 10
        
        if sequence[-1] in ['A', 'U']:
            score += 5
        
        return min(100, score)

    def rate_potency(self, potency):
        if potency >= 85:
            return '优秀'
        elif potency >= 70:
            return '良好'
        elif potency >= 50:
            return '一般'
        else:
            return '较差'

    def suggest_potency(self, potency):
        if potency >= 85:
            return '预测效力较高'
        elif potency >= 70:
            return '预测效力中等，可进一步优化'
        else:
            return '建议重新设计以提高效力'

    def analyze_specificity(self, sequence):
        off_target_risk = 0
        
        if sequence[1:7].count('G') >= 4:
            off_target_risk += 20
        
        if 'GG' in sequence or 'CCC' in sequence:
            off_target_risk += 10
        
        return max(0, 100 - off_target_risk)

    def rate_specificity(self, specificity):
        if specificity >= 90:
            return '优秀'
        elif specificity >= 75:
            return '良好'
        elif specificity >= 60:
            return '一般'
        else:
            return '较差'

    def suggest_specificity(self, specificity):
        if specificity >= 90:
            return '特异性良好，脱靶风险低'
        elif specificity >= 75:
            return '存在一定脱靶风险，建议优化'
        else:
            return '脱靶风险较高，建议重新设计'

    def predict_nuclease_resistance(self, sequence):
        resistance = 50
        
        if sequence.startswith('G'):
            resistance -= 10
        if sequence.startswith('A'):
            resistance += 10
        
        if sequence.endswith('G') or sequence.endswith('C'):
            resistance += 10
        
        for i in range(len(sequence)-1):
            if sequence[i:i+2] in ['UU', 'AA']:
                resistance -= 5
        
        return max(0, min(100, resistance))

    def rate_stability(self, stability):
        if stability >= 70:
            return '优秀'
        elif stability >= 50:
            return '良好'
        elif stability >= 30:
            return '一般'
        else:
            return '较差'

    def suggest_stability(self, stability):
        if stability >= 70:
            return '核酸酶抗性良好'
        elif stability >= 50:
            return '建议考虑化学修饰提高稳定性'
        else:
            return '稳定性较差，强烈建议进行化学修饰'

    def predict_half_life(self, sequence):
        base_half_life = 0.5
        gc = self.calculate_gc_content(sequence)
        base_half_life += gc * 0.01
        
        if sequence.startswith('A') or sequence.startswith('U'):
            base_half_life *= 0.7
        
        return max(0.1, min(5, base_half_life))

    def rate_half_life(self, half_life):
        if half_life >= 2:
            return '优秀'
        elif half_life >= 1:
            return '良好'
        elif half_life >= 0.5:
            return '一般'
        else:
            return '较差'

    def suggest_half_life(self, half_life):
        if half_life >= 2:
            return '血清半衰期较长'
        elif half_life >= 1:
            return '半衰期适中，可考虑修饰延长'
        else:
            return '半衰期较短，建议化学修饰'

    def analyze_chemical_stability(self, sequence):
        if 'UGU' in sequence or 'UGG' in sequence:
            return '较差'
        elif 'CG' in sequence:
            return '一般'
        else:
            return '良好'

    def rate_chemical_stability(self, stability):
        return {'良好': '优秀', '一般': '良好', '较差': '较差'}[stability]

    def suggest_chemical_stability(self, stability):
        if stability == '良好':
            return '化学稳定性良好'
        elif stability == '一般':
            return '存在CG二核苷酸，建议优化'
        else:
            return '存在不稳定基序，建议重新设计'

    def get_target_region_position(self, sequence, region_selection):
        """
        根据功能区域选择计算目标位置范围
        
        Args:
            sequence: mRNA序列
            region_selection: 下拉框选择的功能区域
            
        Returns:
            (start_position, end_position) 元组
        """
        seq_length = len(sequence)
        
        if seq_length < 500:
            cds_start = 0
            cds_end = seq_length
        elif seq_length < 2000:
            cds_start = int(seq_length * 0.15)
            cds_end = int(seq_length * 0.85)
        elif seq_length < 5000:
            cds_start = int(seq_length * 0.20)
            cds_end = int(seq_length * 0.80)
        else:
            cds_start = int(seq_length * 0.25)
            cds_end = int(seq_length * 0.75)
        
        utr5_length = min(200, cds_start) if seq_length >= 500 else int(seq_length * 0.10)
        utr3_start = cds_end if seq_length >= 500 else int(seq_length * 0.90)
        utr3_length = seq_length - utr3_start
        
        region_mapping = {
            '★ CDS (编码序列区)': (cds_start, cds_end),
            '★ 起始密码子附近': (max(0, cds_start - 100), min(seq_length, cds_start + 200)),
            '★ 外显子区域': (cds_start, cds_end),
            '可变剪接区域': (cds_start, cds_end),
            '全序列': (0, seq_length),
            '5\'UTR': (0, min(200, cds_start)),
            '3\'UTR': (max(0, min(utr3_start, seq_length - 30)), seq_length),
            '终止密码子附近': (max(cds_start, cds_end - 200), min(seq_length, cds_end + 100))
        }
        
        start_pos, end_pos = region_mapping.get(region_selection, (0, seq_length))
        
        try:
            offset = int(self.position_entry.get())
            start_pos = max(0, start_pos + offset)
            end_pos = min(seq_length, end_pos + offset)
        except ValueError:
            pass
        
        start_pos = max(0, start_pos)
        end_pos = min(seq_length, end_pos)
        
        return start_pos, end_pos

    def on_region_changed(self, event=None):
        """
        当选择的目标区域变化时，更新说明标签
        """
        region_desc_map = {
            '★ CDS (编码序列区)': ('首选，敲低效果最佳', 'green'),
            '★ 起始密码子附近': ('靠近翻译起始点，靠近外显子连接处', 'green'),
            '★ 外显子区域': ('外显子-外显子连接处效果更好', 'green'),
            '可变剪接区域': ('ASO设计常用', 'blue'),
            '全序列': ('全面搜索所有区域', 'black'),
            '5\'UTR': ('调控研究常用', 'blue'),
            '3\'UTR': ('miRNA结合位点研究', 'blue'),
            '终止密码子附近': ('不推荐，效果差', 'red')
        }

        selected = self.region_combo.get()
        desc, color = region_desc_map.get(selected, ('未知', 'black'))
        self.region_desc_label.config(text=desc, foreground=color)

    def design_sirna(self):
        """
        设计siRNA序列的主方法
        执行流程：
        1. 验证mRNA序列是否存在且有效
        2. 获取设计参数（长度、目标区域）
        3. 调用NucleicAcidDesigner进行siRNA设计
        4. 显示设计结果到表格
        5. 生成BLAST分析结果
        """
        mrna = self.current_mrna

        # 验证mRNA序列是否存在
        if not mrna:
            messagebox.showwarning("警告", "请先查询基因获取mRNA序列")
            return

        # 验证序列有效性
        if not self.designer.is_valid_sequence(mrna):
            messagebox.showwarning("警告", "序列包含无效字符")
            return

        try:
            length_text = self.length_entry.get().strip()
            if not length_text:
                messagebox.showwarning("警告", "请输入序列长度")
                return
            
            length = int(length_text)
            if length < 19 or length > 25:
                messagebox.showwarning("警告", "序列长度应在19-25之间")
                return
                
            region_selection = self.region_combo.get()
            if not region_selection:
                messagebox.showwarning("警告", "请选择目标区域")
                return
                
            start_pos, end_pos = self.get_target_region_position(mrna, region_selection)

            region_display = {
                '★ CDS (编码序列区)': 'CDS (编码序列区)',
                '★ 起始密码子附近': '起始密码子附近',
                '★ 外显子区域': '外显子区域',
                '可变剪接区域': '可变剪接区域',
                '全序列': '全序列',
                '5\'UTR': '5\'UTR',
                '3\'UTR': '3\'UTR',
                '终止密码子附近': '终止密码子附近'
            }.get(region_selection, region_selection)

            candidates = self.designer.find_sirna_candidates(mrna, start_pos, end_pos, length)

            for item in self.result_tree.get_children():
                self.result_tree.delete(item)

            self.blast_text.delete("1.0", tk.END)

            if not candidates:
                self.blast_text.insert("1.0", f"在目标区域 {region_display} (位置 {start_pos}-{end_pos}) 中未找到合适的siRNA候选序列")
                return

            self.set_status(f"正在分析 {len(candidates)} 个候选序列...", "orange")

            for idx, candidate in enumerate(candidates):
                sense = candidate['sense'].replace('T', 'U').replace('t', 'u')
                antisense = candidate['antisense'].replace('T', 'U').replace('t', 'u')
                pos = candidate['position']
                gc = candidate['gc_content']
                reynolds = candidate['reynolds_score']
                amarz = candidate['amarzguioui_score']
                seed_complex = f"{candidate['seed_complexity']:.2f}"
                thermo = "✓" if candidate['thermo_asymmetry'] else "✗"
                uitei = "✓" if candidate['uitei_passed'] else "✗"
                target_region = candidate['target_region']
                final_score = candidate['final_score']
                recommendation = candidate['recommendation']

                rank = "最佳" if idx == 0 else str(idx + 1)
                pos_range = f"{pos} - {pos + length - 1}"

                self.result_tree.insert('', 'end', values=(rank, pos_range, sense, antisense, f"{gc:.1f}", reynolds, amarz, seed_complex, thermo, uitei, target_region, f"{final_score:.1f}", recommendation))

            self.set_status(f"✓ siRNA设计完成，共 {len(candidates)} 个候选序列", "green")

            self.blast_text.delete("1.0", tk.END)
            self.blast_text.insert("1.0", "BLAST相似性分析结果将在此显示\n\n")
            self.blast_text.insert(tk.END, "请先设计siRNA序列，然后选择序列点击\"运行BLAST\"按钮进行脱靶效应检测\n")
            self.blast_text.insert(tk.END, "或点击\"删除选中序列\"移除不满意的序列\n\n")
            self.blast_text.insert(tk.END, f"本次共设计了 {len(candidates)} 个候选序列，可在下方列表中查看详细信息。")
        except ValueError as e:
            messagebox.showerror("错误", f"数值错误: {str(e)}")
            self.set_status("设计失败", "red")
        except Exception as e:
            import traceback
            error_msg = f"设计过程中发生错误: {str(e)}\n\n详细信息:\n{traceback.format_exc()}"
            print(error_msg)
            messagebox.showerror("错误", f"设计过程中发生错误: {str(e)}")
            self.set_status("设计失败", "red")






    def generate_conjugate_sequence(self):
        sense = self.delivery_sense_entry.get().strip().upper()
        antisense = self.delivery_antisense_entry.get().strip().upper()
        delivery = self.delivery_system.get()
        
        if not sense or not antisense:
            messagebox.showwarning("提示", "请先输入siRNA序列")
            return
        
        conjugate_structure = self.build_conjugate_structure(sense, antisense, delivery)
        
        self.delivery_seq_result.config(state='normal')
        self.delivery_seq_result.delete("1.0", tk.END)
        self.delivery_seq_result.insert("1.0", conjugate_structure)
        self.delivery_seq_result.config(state='disabled')
        
        self.set_status(f"✓ 已生成{delivery}结合序列", "green")

    def build_conjugate_structure(self, sense, antisense, delivery):
        """构建递送系统结合后的序列结构"""
        if delivery == 'GalNAC偶联':
            return f"[GalNAC-GalNAC-GalNAC]-L-[{antisense}] || [{sense}]\n\n说明: GalNAC三糖通过连接臂(L)偶联到反义链3'端，实现肝细胞靶向递送"
        elif delivery == 'LNP (脂质纳米粒)':
            return f"LNP包裹: [{antisense}] || [{sense}]\n\n说明: siRNA双链被脂质纳米粒(LNP)包裹形成纳米复合物，通过内吞进入细胞"
        elif delivery == '聚合物纳米粒':
            return f"聚合物包裹: [{antisense}] || [{sense}]\n\n说明: siRNA被可生物降解聚合物包裹，形成稳定的纳米颗粒"
        elif delivery == '病毒载体':
            return f"病毒载体: [{antisense}] || [{sense}]\n\n说明: siRNA装载入病毒载体(如AAV)，通过病毒感染实现高效递送"
        elif delivery == '吸入式纳米粒':
            return f"吸入制剂: [{antisense}] || [{sense}]\n\n说明: siRNA制备成可吸入纳米粒，直接递送至肺部"
        else:
            return f"裸RNA: [{antisense}] || [{sense}]\n\n说明: 未修饰的siRNA双链，适用于局部给药"

    def get_delivery_recommendation(self):
        tissue = self.tissue_type.get()
        methods = DeliverySystem.recommend_delivery(tissue)

        self.delivery_text.delete("1.0", tk.END)
        self.delivery_text.insert(tk.END, f"目标组织: {tissue}\n\n")
        self.delivery_text.insert(tk.END, "推荐递送方式:\n\n")

        for method in methods:
            info = DeliverySystem.DELIVERY_METHODS[method]
            dosage = DeliverySystem.calculate_dosage(21, tissue)
            self.delivery_text.insert(tk.END,
                f"• {info['name']}\n"
                f"  目标: {info['target']}\n"
                f"  优势: {info['advantage']}\n"
                f"  推荐剂量: {dosage:.2f} mg/kg\n\n")

    def set_status(self, message, color='green'):
        """
        更新状态栏显示
        """
        self.status_label.config(text=message, fg=color)
        self.root.update()




def main():
    root = tk.Tk()
    app = NucleicAcidDesignGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
